#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import glob
import argparse

import re
import numpy as np
import matplotlib.pyplot as plt

FILENAME = re.compile('output/.*/pdb/sample_(.*)_(.*)K_.*.pdb')
ATOM = re.compile('ATOM *\d+ *CA ([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')
POTENTIAL = re.compile('REMARK potential: (.*)')


class Residue(object):
    MASS = {
        "ALA ": 71.0779,
        "ARG ": 156.1857,
        "ASN ": 114.1026,
        "ASP ": 115.0874,
        "CYS ": 103.1429,
        "GLU ": 129.114,
        "GLN ": 128.1292,
        "GLY ": 57.0513,
        "HIS ": 137.1393,
        "ILE ": 113.1576,
        "LEU ": 113.1576,
        "LYS ": 128.1723,
        "MET ": 131.1961,
        "PHE ": 147.1739,
        "PRO ": 97.1152,
        "SER ": 87.0773,
        "THR ": 101.1039,
        "TRP ": 186.2099,
        "TYR ": 163.1733,
        "VAL ": 99.1311,
    }

    def __init__(self, amino_acid, index, x, y, z):
        self.amino_acid = amino_acid
        self.index = index
        self.position = np.array((x, y, z))
        self.mass = self.MASS[self.amino_acid]

    def __str__(self):
        return "\t%d: %s %s" % (self.index, self.amino_acid, self.position)


class Protein(object):
    def __init__(self, chain):
        self.residues = []
        self.chain = chain
        self._length = None
        self._radius = None
        self._mass = None
        self._centre_mass = None

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
           # TODO
           # mass
           # centre of mass
           # radius
           # use numpy

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


class Simulation(object):
    def __init__(self, name, chains, conformations):
        self.name = name
        self.conformations = conformations
        self.chains = chains
        self.chains_reversed = {v:k for k,v in self.chains.iteritems()}

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

    def plot_length(self, chain, bins, no_display=False, save=False):
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

        if not no_display:
            plt.show()

        if save:
            # TODO: printable name
            savefile = "%s_%s_length.png" % (self.name.replace('/',''), self.chains[chain])
            plt.savefig(savefile)

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
    s = Simulation.from_glob("UIM/Ub", *args.files, **{"A":"ubiquitin", "B":"UIM"})

    for plot in args.plots:
        plot_method = getattr(s, "plot_%s" % plot, None)
        if not plot_method:
            print "Unsupported plot type: %s" % plot

        for chain in args.chains:
            plot_method(chain, args.bins, args.no_display, args.save)

