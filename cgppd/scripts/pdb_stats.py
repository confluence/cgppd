#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import glob
import os

import re
import math
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import logging

FILENAME = re.compile('output/(.*)/pdb/sample_(.*)_(.*)K_.*.pdb')
SUMMARYFILENAME = re.compile('output/(.*)/summary')
ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')
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


def process_plot(plt, chain_name, plotted_value, args):
    if not args.no_display:
        plt.show()

    if args.save:
        parts = [chain_name, plotted_value]
        if args.log_log:
            parts.append("ll")
        save_name = "_".join(str(p) for p in parts) + ".png"
        save_name = save_name.replace('/', '_').replace(' ', '_')
        plt.savefig(save_name)


def plot_histogram(func):
    def _plot_method(self, chain, args):
        values, graph_name, chain_name, xlabel, ylabel = func(self, chain)
        plotted_value = func.__name__[5:]

        title = "Distribution of %s %s" % (chain_name, graph_name)

        plt.hist(values, bins=args.bins)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        process_plot(plt, chain_name, plotted_value, args)

        plt.close()

    return _plot_method


def plot_vs_sample(func):
    def _plot_method(self, chain, args):
        values, graph_name, chain_name, xlabel, ylabel = func(self, chain)
        plotted_value = func.__name__[5:]

        title = "%s %s vs sample time" % (chain_name, graph_name)

        plt.plot(values, 'bo')
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        process_plot(plt, chain_name, plotted_value, args)

        plt.close()

    return _plot_method


def plot_vs_2powN(func):
    def _plot_method(self, chain, args):
        values, graph_name, chain_name, xlabel, ylabel = func(self, chain)
        plotted_value = func.__name__[5:]

        title = " Root-mean-square %s %s" % (chain_name, graph_name)

        logging.info("Plotting values: %r" % values)

        xvalues = [2**i for i in range(2, len(values) + 2)]

        plt.plot(xvalues, values, 'bo')

        if args.exponent:
            power = reduce(lambda x, y: x/y, [float(n) for n in args.exponent.split('/')])

            fitfunc = lambda p, N: p[0] * np.array(N) ** power # Target function
            errfunc = lambda p, N, y: fitfunc(p, N) - y # Distance to the target function
            p0 = [args.guess_parameter] # Initial guess for the parameters
            p1, success = optimize.leastsq(errfunc, p0[:], args=(xvalues, values))

            print "A = %g" % p1[0]

            plt.plot(xvalues, [fitfunc(p1, x) for x in xvalues], 'b-')

        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if args.log_log:
            plt.xscale('log', basex=2)
            plt.yscale('log')

        process_plot(plt, chain_name, plotted_value, args)

        plt.close()

    return _plot_method


class Simulation(object):
    def __init__(self, name, chains, temperatures):
        self.name = name
        self.chains = chains
        self.temperatures = temperatures
        self._rms_length = {}
        self._rms_radius = {}

    @classmethod
    def write_summary(cls, filenames, summaryfilename, args):
        with open(summaryfilename, 'w') as summaryfile:
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

                summaryfile.write(str(conformation))

    @classmethod
    def from_summary(cls, filename, args):
        with open(filename) as summaryfile:
            summary = summaryfile.read()
            temperatures = defaultdict(list)

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

                temperatures[float(temperature)].append(Conformation(int(sample), float(temperature), float(potential), proteins))

            # TODO: get chains and name from PDB
            chains = dict(n.split(":") for n in args.names)

            name = SUMMARYFILENAME.match(filename).group(1)

            return cls(name, chains, temperatures)

    @property
    def conformations(self):
        return [c for t in sorted(self.temperatures.keys()) for c in self.temperatures[t]]

    @plot_histogram
    def plot_hist_length(self, chain):
        values = [c.proteins[chain].length for c in self.conformations]
        graph_name = "end-to-end length"
        chain_name = "%s %s" % (self.name, self.chains[chain])
        xlabel = u"Length (Å)"
        ylabel = "Frequency (count)"

        return values, graph_name, chain_name, xlabel, ylabel

    @plot_histogram
    def plot_hist_radius(self, chain):
        values = [c.proteins[chain].radius for c in self.conformations]
        graph_name = "radius of gyration"
        chain_name = "%s %s" % (self.name, self.chains[chain])
        xlabel = u"Radius (Å)"
        ylabel = "Frequency (count)"

        return values, graph_name, chain_name, xlabel, ylabel

    @plot_vs_sample
    def plot_sample_length(self, chain):
        values = [c.proteins[chain].length for c in self.conformations]
        graph_name = "end-to-end length"
        chain_name = "%s %s" % (self.name, self.chains[chain])
        xlabel = "Sample"
        ylabel = u"Length (Å)"

        return values, graph_name, chain_name, xlabel, ylabel

    @plot_vs_sample
    def plot_sample_radius(self, chain):
        values = [c.proteins[chain].radius for c in self.conformations]
        graph_name = "radius of gyration"
        chain_name = "%s %s" % (self.name, self.chains[chain])
        xlabel = "Sample"
        ylabel = u"Radius (Å)"

        return values, graph_name, chain_name, xlabel, ylabel

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
        self.name = name
        self.simulations = simulations

    @classmethod
    def from_dirlist(cls, args):
        simulations = []

        for d in args.dirs:
            logging.info("Processing directory '%s'..." % d)

            summaryfilename = "%s/%s/summary" % (args.root_dir, d)
            sim_dir = "%s/%s/pdb" % (args.root_dir, d)

            if not os.path.isfile(summaryfilename):
                if os.path.exists(sim_dir):
                    files = glob.glob("%s/*" % sim_dir)
                    Simulation.write_summary(files, summaryfilename, args)

            if os.path.isfile(summaryfilename):
                s = Simulation.from_summary(summaryfilename, args)
                simulations.append(s)

            if not simulations:
                logging.error("No simulations found. Exiting.")
                sys.exit(1)

        # TODO: fix this
        return cls(args.title, simulations)

    def plot_hist_radius(self, chain, args):
        for s in self.simulations:
            s.plot_hist_radius(chain, args)

    def plot_hist_length(self, chain, args):
        for s in self.simulations:
            s.plot_hist_length(chain, args)

    def plot_sample_radius(self, chain, args):
        for s in self.simulations:
            s.plot_sample_radius(chain, args)

    def plot_sample_length(self, chain, args):
        for s in self.simulations:
            s.plot_sample_length(chain, args)

    @plot_vs_2powN
    def plot_radius(self, chain):
        values = [s.rms_radius[chain] for s in self.simulations]
        graph_name = "radius of gyration"
        chain_name = self.simulations[0].chains[chain]
        xlabel = "Number of residues (count)"
        ylabel = u"Average radius (Å)"

        return values, graph_name, chain_name, xlabel, ylabel

    @plot_vs_2powN
    def plot_length(self, chain):
        values = [s.rms_length[chain] for s in self.simulations]
        graph_name = "end-to-end length"
        chain_name = self.simulations[0].chains[chain]
        xlabel = "Number of residues (count)"
        ylabel = u"Average length (Å)"

        return values, graph_name, chain_name, xlabel, ylabel


AVAILABLE_PLOTS = tuple(n[5:] for n in SimulationSet.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Directories to process", nargs="+")

    parser.add_argument("-r", "--root-dir", default=".", help="Root output directory")
    parser.add_argument("-c", "--chain", action="append", dest="chains", default=[], help="Select a protein by chain label")

    parser.add_argument("-p", "--plot", action="append", dest="plots", default=[], help="Requested plots (available: %s)" % ", ".join(AVAILABLE_PLOTS))
    parser.add_argument("-b", "--bins", help="Number of bins in the histograms", default=30, type=int)
    parser.add_argument("-s", "--save", help="Save plots", action="store_true")
    parser.add_argument("-d", "--no-display", help="Don't display plots", action="store_true")

    parser.add_argument("-v", "--verbose", help="Turn on verbose output", action="store_true")

    # TODO just put these in the files, seriously
    parser.add_argument("-t", "--title", default="Untitled", help="Simulation title")
    parser.add_argument("-n", "--name", action="append", dest="names", default=[], help="Chain name")

    parser.add_argument("-e", "--exponent", help="Power of N for expected scale markers")
    parser.add_argument("-g", "--guess-parameter", type=float, default=1.0, help="Guess constant factor for scale function")
    parser.add_argument("-l", "--log-log", help="Use log-log plot", action="store_true")
    parser.add_argument("-f", "--filter", type=float, help="Filter by temperature")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    if not args.plots:
        logging.error("No plots requested. Exiting.")
        sys.exit(1)

    if not args.chains:
        logging.error("You must select at least one chain. Exiting.")
        sys.exit(1)

    # TODO: add simulation description to config file and write it to PDB file (TITLE)
    # TODO: write chain names into PDB files (COMPND)
    #s = Simulation.from_glob("UIM/Ub", *args.files, **{"A":"ubiquitin", "B":"UIM"})

    s = SimulationSet.from_dirlist(args)

    if s:
        for plot in args.plots:
            plot_method = getattr(s, "plot_%s" % plot, None)
            if not plot_method:
                print "Unsupported plot type: %s" % plot

            for chain in args.chains:
                plot_method(chain, args)

