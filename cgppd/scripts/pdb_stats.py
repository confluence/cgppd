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
    def __init__(self):
        self.proteins = {}

    def add_protein(self, protein):
        self.proteins[protein.chain] = protein

    def __str__(self):
        return "\n\n".join(str(p) for p in self.proteins.values())


class Simulation(object):
    def __init__(self, conformations=None):
        self.conformations = conformations if conformations else []

    @classmethod
    def from_glob(cls, *globs):
        filenames = []
        conformations = []

        for g in globs:
            filenames.extend(glob.glob(g))

        if not filenames:
            print "No files found."
            return

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
                        conformation.add_protein(protein)
                        protein = None

            conformations.append(conformation)

        return cls(conformations)

    def plot_length(self, chain, bins, **kwargs):
        if not all(chain in c.proteins for c in self.conformations):
            return "Chain '%s' does not appear in all conformations." % chain

        lengths = [c.proteins[chain].length for c in self.conformations]

        plt.hist(lengths, bins=bins)
        for k, v in kwargs.iteritems():
            getattr(plt, k)(v)

        plt.show()

    def __str__(self):
        return "\n\n\n".join(str(c) for c in self.simulations[simulation_name])


class SimulationBrowser(object):
    def __init__(self):
        self.simulations = {}

    def add_simulation(self, name, simulation):
        self.simulations[name] = simulation

    def list_simulations(self):
        print "\n".join(self.simulations.keys())


if __name__ == "__main__":
    # TODO: argparse; get list of commands.
    if len(sys.argv) < 2:
        print "No filenames given."
        sys.exit(0)

    sb = SimulationBrowser()
    sb.add_simulation("test", Simulation.from_glob(*sys.argv[1:]))

    ## End-to-end distance

    sb.simulations["test"].plot_length("A", bins=30,
        title="Distribution of ubiquitin end-to-end length in UIM/Ub complex",
        xlabel=u"Length (Å)",
        ylabel="Frequency (count)"
    )

