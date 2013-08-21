#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import glob
import argparse
import time

import re
import math
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
        self._radius = None

    def append_residue(self, residue):
        self.residues.append(residue)

    @property
    def length(self):
        if not self._length:
            start = self.residues[0].position
            end = self.residues[-1].position
            self._length = np.linalg.norm(end - start)

        return self._length

    @property
    def radius(self):
        if not self._radius:
            mass = len(self.residues)
            centre_of_mass = np.array([r.position for r in self.residues]).mean(axis=0)
            # Avoid having to square unnecessary square roots by calculating squared distances with a dot product
            diff_vectors = [(r.position - centre_of_mass) for r in self.residues]
            self._radius = math.sqrt(sum(np.dot(v.T, v) for v in diff_vectors) / mass)

        return self._radius

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


def plot_histogram(func):
    def _plot_method(self, chain, bins, no_display=False, save=False):
        if not all(chain in c.proteins for c in self.conformations):
            raise ValueError("Chain '%s' does not appear in all conformations.", chain)

        values, full_description, xlabel, ylabel = func(self, chain)
        description = func.__name__[5:]

        title = "Distribution of %s %s" % (self.chains[chain], full_description)

        plt.hist(values, bins=bins)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if not no_display:
            plt.show()

        if save:
            # TODO: printable name
            plt.savefig("%s_%s_%s_%d.png" % (self.name.replace('/',''), self.chains[chain], description, time.time()))

        plt.close()

    return _plot_method


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

    @plot_histogram
    def plot_length(self, chain):
        values = [c.proteins[chain].length for c in self.conformations]
        full_description = "end-to-end length"
        xlabel = u"Length (Å)"
        ylabel = "Frequency (count)"

        return values, full_description, xlabel, ylabel


    @plot_histogram
    def plot_radius(self, chain):
        values = [c.proteins[chain].radius for c in self.conformations]
        full_description = "radius of gyration"
        xlabel = u"Radius (Å)"
        ylabel = "Frequency (count)"

        return values, full_description, xlabel, ylabel

    def __str__(self):
        return "\n\n\n".join(str(c) for c in self.simulations[simulation_name])


AVAILABLE_PLOTS = tuple(n[5:] for n in Simulation.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    # TODO: reuse same classes to check potential
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("files", help="Files to process", nargs="+")
    parser.add_argument("-c", "--chain", action="append", dest="chains", default=[], help="Select a protein by chain label")

    parser.add_argument("-p", "--plot", action="append", dest="plots", default=[], help="Requested plots (available: %s)" % ", ".join(AVAILABLE_PLOTS))
    parser.add_argument("-b", "--bins", help="Number of bins in the histograms", default=30, type=int)
    parser.add_argument("-s", "--save", help="Save plots", action="store_true")
    parser.add_argument("-n", "--no-display", help="Don't display plots", action="store_true")

    #parser.add_argument("-v", "--verbose", help="Turn on verbose output", action="store_true")
    args = parser.parse_args()

    if not args.plots:
        # TODO: proper basic logging
        print "No plots requested. Exiting."
        sys.exit(0)

    # TODO: add simulation description to config file and write it to PDB file (TITLE)
    # TODO: write chain names into PDB files (COMPND)
    #s = Simulation.from_glob("UIM/Ub", *args.files, **{"A":"ubiquitin", "B":"UIM"})
    s = Simulation.from_glob("Polyalanine", *args.files, **{"A":"ala9"})

    for plot in args.plots:
        plot_method = getattr(s, "plot_%s" % plot, None)
        if not plot_method:
            print "Unsupported plot type: %s" % plot

        for chain in args.chains:
            plot_method(chain, args.bins, args.no_display, args.save)

