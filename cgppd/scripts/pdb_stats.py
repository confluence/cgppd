#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import time
import glob
import os

import re
import math
import numpy as np
import matplotlib.pyplot as plt

import logging

FILENAME = re.compile('output/(.*)/pdb/sample_(.*)_(.*)K_.*.pdb')
ATOM = re.compile('ATOM *\d+ *CA ([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')
POTENTIAL = re.compile('REMARK potential: (.*)')

RMS_LENGTH_RADIUS = re.compile('name ([^\n]+)\nlength ([A-Z]+ \d+\.\d+)\nradius ([A-Z]+ \d+\.\d+)\n(.*)', re.MULTILINE | re.DOTALL)
FLOAT_PER_CHAIN = re.compile('([A-Z]*) (\d+\.\d+)')
SAMPLE = re.compile('(\d+) at (\d+\.\d+)K with potential (-?\d+\.\d+)\n(.*)', re.MULTILINE | re.DOTALL)
CHAIN = re.compile('chain (.*) length (\d+\.\d+) radius (\d+\.\d+)')

class Residue(object):
    def __init__(self, amino_acid, index, x, y, z):
        self.amino_acid = amino_acid
        self.index = index
        self.position = np.array((x, y, z))

    def __str__(self):
        return "\t%d: %s %s" % (self.index, self.amino_acid, self.position)


class Protein(object):
    def __init__(self, chain, length=None, radius=None):
        self.residues = []
        self.chain = chain
        self._length = length
        self._radius = radius

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
        return "chain %s length %f radius %f\n" % (self.chain, self.length, self.radius)


class Conformation(object):
    def __init__(self, sample, temperature, potential=None, proteins=None):
        self.proteins = proteins or {}
        self.sample = sample
        self.temperature = temperature
        self.potential = potential

    def add_protein(self, protein):
        self.proteins[protein.chain] = protein

    def __str__(self):
        return "sample %d at %fK with potential %f\n%s" % (self.sample, self.temperature, self.potential, "\n".join(str(p) for p in self.proteins.values()))


def plot_histogram(func):
    def _plot_method(self, chain, **kwargs):
        if not all(chain in c.proteins for c in self.conformations):
            raise ValueError("Chain '%s' does not appear in all conformations." % chain)

        values, full_description, xlabel, ylabel = func(self, chain)
        description = func.__name__[5:]

        title = "Distribution of %s %s" % (self.chains[chain], full_description)

        plt.hist(values, bins=kwargs['bins'])
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if not kwargs['no_display']:
            plt.show()

        if kwargs['save']:
            # TODO: printable name
            plt.savefig("hist_%s_%s_%s_%d.png" % (self.name.replace('/',''), self.chains[chain], description, time.time()))

        plt.close()

    return _plot_method


def plot_vs_N(func):
    def _plot_method(self, chain, **kwargs):
        if not all(chain in c.proteins for s in self.simulations for c in s.conformations):
            raise ValueError("Chain '%s' does not appear in all conformations." % chain)

        values, full_description, xlabel, ylabel = func(self, chain)
        description = func.__name__[5:]

        title = " Root-mean-square %s %s" % (self.simulations[0].chains[chain], full_description)

        logging.info("Plotting values: %r" % values)

        plt.plot([i + 2 for i in range(len(values))], values, 'bo')
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        # TODO: how to get this to work
        #plt.set_xscale('log', basex=2)

        if not kwargs['no_display']:
            plt.show()

        if kwargs['save']:
            # TODO: printable name
            plt.savefig("rms_%s_%s_%s_%d.png" % (self.name.replace('/',''), self.simulations[0].chains[chain], description, time.time()))

        plt.close()

    return _plot_method


class Simulation(object):
    def __init__(self, name, chains, conformations, length=None, radius=None):
        self.name = name
        self.conformations = conformations
        self.chains = chains
        self._rms_length = length or {}
        self._rms_radius = radius or {}

    @classmethod
    def from_filelist(cls, filenames):
        conformations = []

        for filename in filenames:
            with open(filename) as pdbfile:
                logging.debug("Parsing file '%s'..." % filename)
                name, sample, temperature = FILENAME.match(filename).groups()
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

        # TODO: get chains and name from PDB
        chains = {"A": name}

        return cls(name, chains, conformations)

    @classmethod
    def from_summary(cls, filename):
        with open(filename) as summaryfile:
            summary = summaryfile.read()

            name, length, radius, summary = RMS_LENGTH_RADIUS.search(summary).groups()

            length_dict = dict((k, float(v)) for (k, v) in FLOAT_PER_CHAIN.findall(length))
            radius_dict = dict((k, float(v)) for (k, v) in FLOAT_PER_CHAIN.findall(radius))

            conformations = []

            for s in re.split("sample", summary):
                if not s:
                    continue

                sample, temperature, potential, chains = SAMPLE.search(s).groups()
                proteins = {}

                for c in chains.split("\n"):
                    if not c:
                        continue
                    chain, length, radius = CHAIN.match(c).groups()
                    proteins[chain] = Protein(chain, float(length), float(radius))

                conformations.append(Conformation(int(sample), float(temperature), float(potential), proteins))

            # TODO: get chains and name from PDB
            chains = {"A": name}

            return cls(name, chains, conformations, length_dict, radius_dict)


    def write_summary(self, filename):
        with open(filename, 'w') as summaryfile:
            summaryfile.write(str(self))

            for c in self.conformations:
                summaryfile.write(str(c))

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
        length = " ".join("%s %f" % (k, v) for (k, v) in self.rms_length.iteritems())
        radius = " ".join("%s %f" % (k, v) for (k, v) in self.rms_radius.iteritems())
        return "name %s\nlength %s\nradius %s\n" % (self.name, length, radius)

    @property
    def rms_radius(self):
        if not self._rms_radius:
            for chain in self.chains:
                values = np.array([c.proteins[chain].radius for c in self.conformations])
                self._rms_radius[chain] = np.sqrt(np.mean(values**2))

        return self._rms_radius

    @property
    def rms_length(self):
        if not self._rms_length:
            for chain in self.chains:
                values = np.array([c.proteins[chain].length for c in self.conformations])
                self._rms_length[chain] = np.sqrt(np.mean(values**2))

        return self._rms_length


class SimulationSet(object):
    def __init__(self, name, simulations):
        self.simulations = simulations

    @classmethod
    def from_dirlist(cls, name, *dirs):
        simulations = []

        for d in dirs:
            logging.info("Processing directory '%s'..." % d)

            summaryfilename = "%s/summary" % d

            if os.path.isfile(summaryfilename):
                s = Simulation.from_summary(summaryfilename)

            else:
                files = glob.glob("%s/pdb/*" % d)
                s = Simulation.from_filelist(files)
                s.write_summary(summaryfilename)

            simulations.append(s)

        return cls(name, simulations)

    def plot_hist_radius(self, chain, **kwargs):
        for s in self.simulations:
            s.plot_radius(chain, **kwargs)

    def plot_hist_length(self, chain, **kwargs):
        for s in self.simulations:
            s.plot_length(chain, **kwargs)

    @plot_vs_N
    def plot_radius(self, chain):
        values = [s.rms_radius[chain] for s in self.simulations]
        full_description = "radius of gyration"
        xlabel = u"Average radius (Å)"
        ylabel = "Number of residues (count)"

        return values, full_description, xlabel, ylabel

    @plot_vs_N
    def plot_length(self, chain):
        values = [s.rms_length[chain] for s in self.simulations]
        full_description = "end-to-end length"
        xlabel = u"Average length (Å)"
        ylabel = "Number of residues (count)"

        return values, full_description, xlabel, ylabel


AVAILABLE_PLOTS = tuple(n[5:] for n in SimulationSet.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Directories to process", nargs="+")
    parser.add_argument("-c", "--chain", action="append", dest="chains", default=[], help="Select a protein by chain label")

    parser.add_argument("-p", "--plot", action="append", dest="plots", default=[], help="Requested plots (available: %s)" % ", ".join(AVAILABLE_PLOTS))
    parser.add_argument("-b", "--bins", help="Number of bins in the histograms", default=30, type=int)
    parser.add_argument("-s", "--save", help="Save plots", action="store_true")
    parser.add_argument("-n", "--no-display", help="Don't display plots", action="store_true")

    parser.add_argument("-v", "--verbose", help="Turn on verbose output", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    if not args.plots:
        # TODO: proper basic logging
        print "No plots requested. Exiting."
        sys.exit(0)

    # TODO: add simulation description to config file and write it to PDB file (TITLE)
    # TODO: write chain names into PDB files (COMPND)
    #s = Simulation.from_glob("UIM/Ub", *args.files, **{"A":"ubiquitin", "B":"UIM"})

    s = SimulationSet.from_dirlist("Polyalanine", *args.dirs)

    if s:
        for plot in args.plots:
            plot_method = getattr(s, "plot_%s" % plot, None)
            if not plot_method:
                print "Unsupported plot type: %s" % plot

            for chain in args.chains:
                plot_method(chain, bins=args.bins, no_display=args.no_display, save=args.save)

