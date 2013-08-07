#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

ATOM = re.compile('ATOM *\d+ *CA ([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')

class Residue(object):
    def __init__(self, amino_acid, index, x, y, z):
        self.amino_acid = amino_acid
        self.index = index
        self.position = np.array((x, y, z))

    def __str__(self):
        return "\t%d: %s %s" % (self.index, self.amino_acid, self.position)


class Protein(object):
    def __init__(self, chain):
        self.residues = []
        self.chain = chain

    def append_residue(self, residue):
        self.residues.append(residue)

    def __str__(self):
        return "%s:\n%s" % (self.chain, "\n".join([str(r) for r in self.residues]))


class Conformation(object):
    def __init__(self):
        self.proteins = []

    def append_protein(self, protein):
        self.proteins.append(protein)

    def __str__(self):
        return "\n\n".join([str(p) for p in self.proteins])

if __name__ == "__main__":
    # TODO: proper argparse
    if len(sys.argv) < 2:
        print "No filenames given."
        sys.exit(0)

    filenames = []

    for arg in sys.argv[1:]:
        filenames.extend(glob.glob(arg))

    if not filenames:
        print "No files found."
        sys.exit(0)

    conformations = []

    for filename in filenames:
        with open(filename) as pdbfile:
            conformation = Conformation()
            protein = None

            for line in pdbfile:
                if line.startswith("ATOM"):
                    amino_acid, chain, index, x, y, z = ATOM.match(line).groups()
                    if not protein:
                        protein = Protein(chain)
                    protein.append_residue(Residue(amino_acid, int(index), float(x), float(y), float(z)))

                elif line.startswith("TER"):
                    conformation.append_protein(protein)
                    protein = None

        conformations.append(conformation)

    #print "\n\n\n".join([str(c) for c in conformations])

    # End-to-end distance

    distances = {}

    for conformation in conformations:
        for protein in conformation.proteins:
            start = protein.residues[0].position
            end = protein.residues[-1].position
            distance = np.linalg.norm(end - start)

            if protein.chain not in distances:
                distances[protein.chain] = []

            distances[protein.chain].append(distance)


    plt.hist(distances["A"], bins=30)
    plt.title("Distribution of ubiquitin end-to-end length in UIM/Ub complex")
    plt.xlabel(u"Length (Å)")
    plt.ylabel("Frequency")
    plt.show()