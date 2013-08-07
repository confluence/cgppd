#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

FILENAME = re.compile('output/.*/pdb/sample_(.*)_(.*)K_.*.pdb')
ATOM = re.compile('ATOM *\d+ *CA ([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')
POTENTIAL = re.compile('REMARK potential: (.*)')

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
        self._length = None

    def append_residue(self, residue):
        self.residues.append(residue)

    @property
    def length(self):
        if not self._length:
            start = self.residues[0].position
            end = self.residues[-1].position
            self._length = np.linalg.norm(end - start)

        return self._length

    def __str__(self):
        return "%s:\n%s" % (self.chain, "\n".join(str(r) for r in self.residues))


class Conformation(object):
    def __init__(self, sample, temperature, potential=None):
        self.proteins = {}
        self.sample = sample
        self.temperature = temperature
        self.potential = potential

    def add_protein(self, protein):
        self.proteins[protein.chain] = protein

    def __str__(self):
        return "Sample %d at %gK with potential %g\n%s" % (self.sample, self.temperature, self.potential, "\n".join(str(p) for p in self.proteins.values()))


class Simulation(object):
    def __init__(self, name, chains, conformations):
        self.name = name
        self.conformations = conformations
        self.chains = chains

    @classmethod
    def from_glob(cls, name, *globs, **chains):
        filenames = []
        conformations = []

        for g in globs:
            filenames.extend(glob.glob(g))

        if not filenames:
            print "No files found."
            return

        for filename in filenames:
            with open(filename) as pdbfile:
                sample, temperature = FILENAME.match(filename).groups()
                conformation = Conformation(int(sample), float(temperature))

                protein = None

                for line in pdbfile:
                    if line.startswith("ATOM"):
                        amino_acid, chain, index, x, y, z = ATOM.match(line).groups()
                        if not protein:
                            protein = Protein(chain)
                        protein.append_residue(Residue(amino_acid, int(index), float(x), float(y), float(z)))

                    elif line.startswith("TER"):
                        conformation.add_protein(protein)
                        protein = None

                    elif line.startswith("REMARK potential"):
                        conformation.potential = float(POTENTIAL.match(line).group(1))

            conformations.append(conformation)

        return cls(name, chains, conformations)

    def plot_length(self, chain, bins, save=False):
        if not all(chain in c.proteins for c in self.conformations):
            return "Chain '%s' does not appear in all conformations." % chain

        lengths = [c.proteins[chain].length for c in self.conformations]
        x_axis = u"Length (Å)"
        y_axis = "Frequency (count)"
        title = "Distribution of %s end-to-end length" % self.chains[chain]

        plt.hist(lengths, bins=bins)
        plt.title(title)
        plt.xlabel(x_axis)
        plt.ylabel(y_axis)

        if save:
            # TODO: printable name
            savefile = "%s_%s_length.png" % (self.name.replace('/',''), self.chains[chain])
            plt.savefig(savefile)
        else:
            plt.show()

    def __str__(self):
        return "\n\n\n".join(str(c) for c in self.simulations[simulation_name])


if __name__ == "__main__":
    # TODO: argparse; get list of commands.
    if len(sys.argv) < 2:
        print "No filenames given."
        sys.exit(0)

    # TODO: add simulation description to config file and write it to PDB file
    # TODO: write chain names into PDB files
    s = Simulation.from_glob("UIM/Ub", *sys.argv[1:], **{"A":"ubiquitin", "B":"UIM"})

    print s.conformations[0]

    ## End-to-end distance

    s.plot_length("A", bins=30)

