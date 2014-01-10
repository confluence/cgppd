#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import glob
import os

import re
import math
import itertools
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
        self.residues.append(residue)

    def measure(self):
        if self.length is None:
            start = self.residues[0].position
            end = self.residues[-1].position
            self.length = np.linalg.norm(end - start)

        if self.radius is None:
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
    FILENAME = re.compile('.*/pdb/sample_(.*)_(.*)K_ *(.*).pdb')
    NAME = re.compile('(\d*)$')
    ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')
    SAMPLE = re.compile('(\d+) at (\d+\.\d+)K with potential (-?\d+\.\d+)\n(.*)', re.MULTILINE | re.DOTALL)
    CHAIN = re.compile('chain (.*) length (\d+\.\d+) radius (\d+\.\d+)')

    def __init__(self, temperatures, name, N=None):
        self.temperatures = temperatures # actually needs to be a dict
        self.name = name
        self.N = N

    @classmethod
    def write_summary(cls, filenames, summaryfilename, args):
        with open(summaryfilename, 'w') as summaryfile:
            for filename in filenames:
                with open(filename) as pdbfile:
                    logging.debug("Parsing file '%s'..." % filename)
                    sample_id, temperature, potential = cls.FILENAME.match(filename).groups()

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
                    for chain in sample.chains:
                        chain.measure()
                        summaryfile.write(str(chain))

    @classmethod
    def from_summary(cls, filename, args):
        name = filename.split('/')[-2].split("_")[0]
        N = cls.NAME.search(name).group(1) or None

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

            return cls(temperatures, name, N)

    @property
    def temperatures_list(self):
        return [self.temperatures[t] for t in sorted(self.temperatures.keys())]


class Dataset(object):
    def __init__(self, simulation_set, args):
        self.data = {}

        for sim in simulation_set.simulations:
            for t in sim.temperatures_list:
                if args.temperatures and t not in args.temperatures:
                    continue

                for s in t.samples:
                    for c in s.chains:
                        if args.chains and c.chain_id not in args.chains:
                            continue

                        key = (sim.N or sim.name, c.chain_id, c.temperature)

                        if key not in self.data:
                            self.data[key] = []
                        self.data[key].append(c)

        self.simulations, self.chain_ids, self.temperatures = [sorted(list(set(l))) for l in zip(*self.data.keys())]

        num_t = len(self.temperatures)
        dimy = np.ceil(num_t**0.5)
        dimx = np.ceil(num_t/dimy)
        self.dimx, self.dimy = int(dimx), int(dimy)

        self.measurements = set(simulation_set.MEASUREMENTS)
        if args.measurements:
            self.measurements &= set(args.measurements)


class SimulationSet(object):
    MEASUREMENTS = ('length', 'radius')

    def __init__(self, simulations):
        self.simulations = simulations # list

    @classmethod
    def from_dirlist(cls, args):
        simulations = []

        for d in args.dirs:
            logging.info("Processing directory '%s'..." % d)

            summaryfilename = "%s/summary" % d
            sim_dir = "%s/pdb" % d

            if not os.path.isfile(summaryfilename):
                if os.path.exists(sim_dir):
                    files = glob.glob("%s/*" % sim_dir)
                    Simulation.write_summary(files, summaryfilename, args)

            if os.path.isfile(summaryfilename):
                s = Simulation.from_summary(summaryfilename, args)
                simulations.append(s)

        return cls(simulations)

    def _display_plot(self, plt, fig, title, xlabel, ylabel):
        plt.tight_layout()
        fig.suptitle(title, fontsize=14)
        fig.text(0.5, 0.02, xlabel, ha='center', va='center')
        fig.text(0.02, 0.5, ylabel, ha='center', va='center', rotation='vertical')
        plt.show()
        plt.close()

    def plot_histogram(self, dataset, args):
        for measurement in dataset.measurements:

            for sim in dataset.simulations:
                for chain_id in dataset.chain_ids:
                    fig = plt.figure()

                    for i, temperature in enumerate(dataset.temperatures):
                        ax = fig.add_subplot(dataset.dimy, dataset.dimx, i + 1)

                        chains = dataset.data[(sim, chain_id, temperature)]

                        values = np.array([getattr(c, measurement) for c in chains])
                        ax.hist(values)

                        ax.set_title("%sK" % temperature)

                        labels = ax.get_xticklabels()
                        for label in labels:
                            label.set_rotation(30)
                            label.set_fontsize(10)

                    title = "%s %s" % (chain_id, measurement)
                    xlabel = u"%s (Å)" % measurement
                    ylabel = "Frequency (count)"


                    self._display_plot(plt, fig, title, xlabel, ylabel)

    def plot_samples(self, dataset, args):
        for measurement in dataset.measurements:

            for sim in dataset.simulations:
                for chain_id in dataset.chain_ids:
                    fig = plt.figure()

                    for i, temperature in enumerate(dataset.temperatures):
                        ax = fig.add_subplot(dataset.dimy, dataset.dimx, i + 1)

                        chains = dataset.data[(sim, chain_id, temperature)]

                        values = np.array([getattr(c, measurement) for c in chains])
                        ax.plot(values, 'bo')

                        ax.set_title("%sK" % temperature)

                        labels = ax.get_xticklabels()
                        for label in labels:
                            label.set_rotation(30)
                            label.set_fontsize(10)

                    title = "%s %s" % (chain_id, measurement)
                    xlabel = "Sample time"
                    ylabel = u"%s (Å)" % measurement

                    self._display_plot(plt, fig, title, xlabel, ylabel)

    def plot_vs_N(self, dataset, args):
        try:
            xvalues = np.array([int(N) for N in dataset.simulations])
        except TypeError:
            logging.error("This dataset does not have a range of residue sizes.")
            sys.exit(1)

        for measurement in dataset.measurements:
            for chain_id in dataset.chain_ids:
                values_per_temperature = []

                if args.aggregate: # aggregate all temperatures
                    values = []

                    for N in dataset.simulations:
                        chains = [c for t in dataset.temperatures for c in dataset.data[(N, chain_id, t)]]
                        chain_values = np.array([getattr(c, measurement) for c in chains])
                        values.append(np.sqrt(np.mean(chain_values**2)))

                    values_per_temperature.append((None, values))

                else:
                    for temperature in dataset.temperatures:
                        values = []

                        for N in dataset.simulations:
                            chains = dataset.data[(N, chain_id, temperature)]
                            chain_values = np.array([getattr(c, measurement) for c in chains])
                            values.append(np.sqrt(np.mean(chain_values**2)))

                        values_per_temperature.append((temperature, values))

                fig = plt.figure()

                for i, (temperature, values) in enumerate(values_per_temperature):
                    if temperature is not None:
                        ax = fig.add_subplot(dataset.dimy, dataset.dimx, i + 1)
                    else:
                        ax = fig.add_subplot(1, 1, i + 1)

                    plt.plot(xvalues, values, 'bo')

                    if args.fit:
                        power = reduce(lambda x, y: x/y, [float(n) for n in args.fit.split('/')])

                        fitfunc = lambda p, N: p[0] * np.array(N) ** power # Target function
                        errfunc = lambda p, N, y: fitfunc(p, N) - y # Distance to the target function
                        p0 = [1.0] # Initial guess for the parameters
                        p1, success = optimize.leastsq(errfunc, p0[:], args=(xvalues, values))

                        print "A = %g" % p1[0]

                        ax.plot(xvalues, [fitfunc(p1, x) for x in xvalues], 'b-')

                        if temperature is not None:
                            ax.set_title("%sK" % temperature)
                        else:
                            ax.set_title("All temperatures (%sK - %sK)" % (dataset.temperatures[0], dataset.temperatures[-1]))

                        ax.set_xscale('log', basex=2) # TODO: investigate using the loglog function instead; maybe add an option for it
                        ax.set_yscale('log')

                title = "%s root-mean-square %s" % (chain_id, measurement)
                xlabel = "Number of residues (count)"
                ylabel = u"Average %s (Å)" % measurement

                self._display_plot(plt, fig, title, xlabel, ylabel)

    def all_plots(self, args):
        dataset = Dataset(self, args)

        for plot in args.plots:
            plot_method = getattr(self, "plot_%s" % plot, None)

            if not plot_method:
                print "Unsupported plot type: %s" % plot
                sys.exit(1)

            plot_method(dataset, args)


PLOTS = tuple(n[5:] for n in SimulationSet.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")

    parser.add_argument("-p", "--plot", action="append", dest="plots", default=[], help="Select plot(s) (available: %s)" % ", ".join(PLOTS))
    parser.add_argument("-m", "--measurement", action="append", dest="measurements", default=[], help="Select measurement(s) (available: %s)" % ", ".join(SimulationSet.MEASUREMENTS))
    parser.add_argument("-c", "--chain", action="append", dest="chains", default=[], help="Select chain(s)")
    parser.add_argument("-t", "--temperature", action="append", type=float, dest="temperatures", default=[], help="Select temperature(s)")
    parser.add_argument("-a", "--aggregate", help="Aggregate temperatures (plot vs N only)", action="store_true")
    parser.add_argument("-f", "--fit", help="Attempt to fit values to function with the provided exponent (plot vs N only)")

    parser.add_argument("-v", "--verbose", help="Turn on verbose output", action="store_true")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    s = SimulationSet.from_dirlist(args)
    s.all_plots(args)



