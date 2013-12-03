#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import glob
import os

import re
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import logging


class Residue(object):
    def __init__(self, amino_acid, index, x, y, z):
        self.amino_acid = amino_acid
        self.index = index
        self.position = np.array((x, y, z))

    def __str__(self):
        return "%d: %s %s" % (self.index, self.amino_acid, self.position)


class Chain(object):
    def __init__(self, chain_id, length=None, radius=None, temperature=None, residues=None):
        self.chain_id = chain_id
        self.length = length
        self.radius = radius
        self.temperature = temperature
        self.residues = residues or []

    def add_residue(self, residue):
        residues.append(residue)

    def measure(self):
        if length is None:
            start = self.residues[0].position
            end = self.residues[-1].position
            self.length = np.linalg.norm(end - start)

        if radius is None:
            mass = len(self.residues)
            centre_of_mass = np.array([r.position for r in self.residues]).mean(axis=0)
            # Avoid having to square unnecessary square roots by calculating squared distances with a dot product
            diff_vectors = [(r.position - centre_of_mass) for r in self.residues]
            self.radius = math.sqrt(sum(np.dot(v.T, v) for v in diff_vectors) / mass)


    def __str__(self):
        return "chain %s length %f radius %f" % (self.chain_id, self.length, self.radius)


class Sample(object):
    def __init__(self, sample_id, temperature, potential, chains=None):
        self.sample_id = sample_id
        self.temperature = temperature
        self.potential = potential
        self.chains = chains or []

    def add_chain(self, chain):
        self.chains.append(chain)

    def __str__(self):
        return "sample %d at %fK with potential %f" % (self.sample_id, self.temperature, self.potential)


class Temperature(object):
    def __init__(self, temperature, samples=None):
        self.temperature = temperature
        self.samples = samples or []

    def add_sample(self, sample):
        self.samples.append(sample)

    def __str__(self):
        return "temperature %f with %d samples" % (self.temperature, len(self.samples))


class Simulation(object):
    FILENAME = re.compile('output/(.*)/pdb/sample_(.*)_(.*)K_ *(.*).pdb')
    ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')
    SAMPLE = re.compile('(\d+) at (\d+\.\d+)K with potential (-?\d+\.\d+)\n(.*)', re.MULTILINE | re.DOTALL)
    CHAIN = re.compile('chain (.*) length (\d+\.\d+) radius (\d+\.\d+)')

    def __init__(self, temperatures):
        self.temperatures = temperatures # actually needs to be a dict

    @classmethod
    def write_summary(cls, filenames, summaryfilename, args):
        with open(summaryfilename, 'w') as summaryfile:
            for filename in filenames:
                with open(filename) as pdbfile:
                    logging.debug("Parsing file '%s'..." % filename)
                    name, sample_id, temperature, potential = cls.FILENAME.match(filename).groups()

                    sample = Sample(int(sample_id), float(temperature), float(potential))
                    chain = None

                    for line in pdbfile:
                        if line.startswith("ATOM"):
                            amino_acid, chain_id, index, x, y, z = cls.ATOM.match(line).groups()

                            if not chain:
                                chain = Chain(chain_id)
                            chain.add_residue(Residue(amino_acid, int(index), float(x), float(y), float(z)))

                        elif line.startswith("TER"):
                            sample.add_chain(chain)
                            chain = None

                    summaryfile.write(str(sample))
                    for chain in sample.chains.itervalues():
                        chain.measure()
                        summaryfile.write(str(chain))

    @classmethod
    def from_summary(cls, filename, args):
        with open(filename) as summaryfile:
            summary = summaryfile.read()
            temperatures = {}

            for s in re.split("sample", summary):
                if not s:
                    continue

                sample_id, temperature, potential, chain_lines = cls.SAMPLE.search(s).groups()
                sample_id, temperature, potential = int(sample_id), float(temperature), float(potential)

                sample = Sample(sample_id, temperature, potential)

                for c in chain_lines.split("\n"):
                    if not c:
                        continue
                    chain_id, length, radius = cls.CHAIN.match(c).groups()
                    sample.add_chain(Chain(chain_id, float(length), float(radius), temperature))

                if temperature not in temperatures:
                    temperatures[temperature] = Temperature(temperature)

                temperatures[temperature].add_sample(sample)

            return cls(temperatures)

    @property
    def temperatures_list(self):
        return [self.temperatures[t] for t in sorted(self.temperatures.keys())]


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


class SimulationSet(object):
    def __init__(self, simulations):
        self.simulations = simulations # list

    @classmethod
    def from_dirlist(cls, args):
        simulations = []

        if args.directory_prefix:
            # Automatically select and order all directories in the root directory
            dirs = natural_sort([d for d in os.listdir(args.root_dir) if d.startswith(args.directory_prefix)])
        else:
            dirs = args.dirs

        for d in dirs:
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

        return cls(simulations)

    def _display_or_save_plot(self, figure):
        if not args.no_display:
            display(figure)

        if args.save:
            figure.savefig("changemyname.png")

    def plot_histogram(self, chains, args):
        # split up by chain and measurement and temperature in here?
        # maybe split it up outside and pass it in?
        pass

    def plot_samples(self, chains, args):
        # split up by chain and measurement and temperature in here
        pass

    def plot_vs_N(self, chains, args):
        # split up by chain and measurement and temperature in here
        # aggregate (or not) in here
        pass

    def all_plots(self, args):
        # filter stuff here
        chains = [c for sim in self.simulations for t in sim.temperatures for s in t.samples for c in s.chains]

        if args.chains:
            chains = [c for c in chains if c.chain_id in args.chains]

        if args.temperatures:
            chains = [c for c in chains if c.temperature in args.temperatures]

        for plot in args.plots:
            plot_method = getattr(self, "plot_%s" % plot, None)

            if not plot_method:
                print "Unsupported plot type: %s" % plot

            plot_method(args)


PLOTS = tuple(n[5:] for n in SimulationSet.__dict__ if n.startswith("plot_"))
MEASUREMENTS = ('length', 'radius')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")

    parser.add_argument("-r", "--root-dir", default=".", help="Root directory")
    parser.add_argument("-d", "--directory-prefix", help="Directory prefix")

    parser.add_argument("-p", "--plot", action="append", dest="plots", default=[], help="Select plot(s) (available: %s)" % ", ".join(PLOTS))
    parser.add_argument("-m", "--measurement", action="append", dest="measurements", default=[], help="Select measurement(s) (available: %s)" % ", ".join(MEASUREMENTS))
    parser.add_argument("-c", "--chain", action="append", dest="chains", default=[], help="Select chain(s)")
    parser.add_argument("-t", "--temperature", action="append", dest="temperatures", default=[], help="Select temperature(s)")
    parser.add_argument("-a", "--aggregate", help="Aggregate temperatures (plot vs N only)", action="store_true")

    parser.add_argument("-v", "--verbose", help="Turn on verbose output", action="store_true")
    parser.add_argument("-s", "--save", help="Save plots", action="store_true")
    parser.add_argument("-n", "--no-display", help="Don't display plots", action="store_true")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    s = SimulationSet.from_dirlist(args)
    s.all_plots(args)



