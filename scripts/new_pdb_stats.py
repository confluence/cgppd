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
        return "chain %s length %f radius %f\n" % (self.chain_id, self.length, self.radius)


class Sample(object):
    def __init__(self, sample_id, temperature, potential, chains=None):
        self.sample_id = sample_id
        self.temperature = temperature
        self.potential = potential
        self.chains = chains or []

    def add_chain(self, chain):
        self.chains.append(chain)

    def __str__(self):
        return "sample %d at %fK with potential %f\n" % (self.sample_id, self.temperature, self.potential)


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
    REMARK_SAMPLE = re.compile('REMARK sample: (\d+)')
    REMARK_POTENTIAL = re.compile('REMARK potential: (-?\d+\.\d+)')
    FRAME = re.compile('TITLE *frame t= -?(\d+).000')
    
    SAMPLE = re.compile('(\d+) at (\d+\.\d+)K with potential (-?\d+\.\d+)\n(.*)', re.MULTILINE | re.DOTALL)
    CHAIN = re.compile('chain (.*) length (\d+\.\d+) radius (\d+\.\d+)')

    def __init__(self, temperatures, name, N=None):
        self.temperatures = temperatures # actually needs to be a dict
        self.name = name
        self.N = N

    @classmethod
    def write_summary_from_trajectory(cls, trajectoryfilename, summaryfilename, args):
        with open(trajectoryfilename, 'r') as trajectoryfile:
            with open(summaryfilename, 'w') as summaryfile:
                logging.debug("Parsing trajectory file '%s' ..." % (trajectoryfilename))

                # we probably need to postprocess with cgppd to find the potentials -- or do we?
                # also probably need to store information about the clusters, e.g. number of samples per cluster
                # we need to measure the length of the whole molecule, not just the chains
                sample = Sample(0, 300.0, 0)
                chain = None
                
                for line in trajectoryfile:
                    if line.startswith("REMARK sample:"): # unclustered trajectory
                        sample.sample_id = int(cls.REMARK_SAMPLE.match(line).group(1))
                    elif line.startswith("REMARK potential:"): # unclustered trajectory
                        sample.potential = float(cls.REMARK_POTENTIAL.match(line).group(1))
                    elif "frame" in line: # clustered trajectory
                        sample.sample_id = int(cls.FRAME.match(line).group(1))
                    elif line.startswith("ATOM"):
                        amino_acid, chain_id, index, x, y, z = cls.ATOM.match(line).groups()

                        if not chain:
                            chain = Chain(chain_id)
                        chain.add_residue(Residue(amino_acid, int(index), float(x), float(y), float(z)))

                    elif line.startswith("TER"): # end of chain
                        sample.add_chain(chain)
                        chain = None
                    
                    elif line.startswith("ENDMDL"): # end of sample
                        summaryfile.write(str(sample))

                        # quick and dirty hack to calculate stats for multi-chain molecule
                        # length relies on sensible ordering of chains in the PDB file
                        molecule = Chain("all")
                        
                        for chain in sample.chains:
                            chain.measure()
                            summaryfile.write(str(chain))

                            molecule.residues.extend(chain.residues)

                        molecule.measure()
                        summaryfile.write(str(molecule))
                            
                        sample = Sample(0, 300.0, 0)
                        chain = None
                            
                

    @classmethod
    def write_summary(cls, filenames, summaryfilename, args):
        num_files = len(filenames)
        with open(summaryfilename, 'w') as summaryfile:
            for i, filename in enumerate(filenames):
                with open(filename) as pdbfile:
                    logging.debug("Parsing file '%s' (%d of %d)..." % (filename, i, num_files))
                    if not i % 1000:
                        logging.info("%d of %d files parsed..." % (i, num_files))
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
        directory, basename = os.path.split(filename)
        if basename == "summary":
            name = os.path.split(directory)[1]
        else:
            name = os.path.splitext(basename)[0].split("_")[1]
        
        name = os.path.basename(filename).split("_")[0]
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
                if args.temperatures and t.temperature not in args.temperatures:
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
    def from_trajectorylist(cls, args):
        simulations = []

        for trajectoryfilename in args.trajectories:
            logging.info("Processing trajectory '%s'..." % trajectoryfilename)

            summaryfilename = "summary_%s" % trajectoryfilename

            if not os.path.isfile(summaryfilename):
                logging.info("Writing new summary file...")
                Simulation.write_summary_from_trajectory(trajectoryfilename, summaryfilename, args)

            if os.path.isfile(summaryfilename):
                logging.info("Loading summary file...")
                s = Simulation.from_summary(summaryfilename, args)
                simulations.append(s)

        return cls(simulations)

    @classmethod
    def from_dirlist(cls, args):
        simulations = []

        for d in args.dirs:
            logging.info("Processing directory '%s'..." % d)

            summaryfilename = "%s/summary" % d
            sim_dir = "%s/pdb" % d

            if not os.path.isfile(summaryfilename):
                logging.info("Writing new summary file...")
                if os.path.exists(sim_dir):
                    files = glob.glob("%s/*" % sim_dir)
                    Simulation.write_summary(files, summaryfilename, args)

            if os.path.isfile(summaryfilename):
                logging.info("Loading summary file...")
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

                for temperature in dataset.temperatures:
                    values = []

                    for N in dataset.simulations:
                        chains = dataset.data[(N, chain_id, temperature)]
                        chain_values = np.array([getattr(c, measurement) for c in chains])
                        values.append(np.sqrt(np.mean(chain_values**2)))

                    values_per_temperature.append((temperature, values))

                fig = plt.figure()

                for i, (temperature, values) in enumerate(values_per_temperature):
                    ax = fig.add_subplot(dataset.dimy, dataset.dimx, i + 1)

                    plt.plot(xvalues, values, 'bo')

                    if args.fit:
                        power = reduce(lambda x, y: x/y, [float(n) for n in args.fit.split('/')])

                        fitfunc = lambda p, N: p[0] * np.array(N) ** power # Target function
                        errfunc = lambda p, N, y: fitfunc(p, N) - y # Distance to the target function
                        p0 = [1.0] # Initial guess for the parameters
                        p1, success = optimize.leastsq(errfunc, p0[:], args=(xvalues, values))

                        logging.debug("A = %g" % p1[0])

                        ax.plot(xvalues, [fitfunc(p1, x) for x in xvalues], 'b-')

                        ax.set_title("%sK" % temperature)
  

                        ax.set_xscale('log', basex=2) # TODO: investigate using the loglog function instead; maybe add an option for it
                        ax.set_yscale('log')

                title = "%s root-mean-square %s" % (chain_id, measurement)
                xlabel = "Number of residues (count)"
                ylabel = u"Average %s (Å)" % measurement

                self._display_plot(plt, fig, title, xlabel, ylabel)

    def all_plots(self, args):
        dataset = Dataset(self, args)

        if not args.plots:
            logging.info("No plots selected. Preprocessing only.") 

        for plot in args.plots:
            plot_method = getattr(self, "plot_%s" % plot, None)

            if not plot_method:
                logging.warn("Unsupported plot type: %s" % plot)
                sys.exit(1)

            plot_method(dataset, args)


PLOTS = tuple(n[5:] for n in SimulationSet.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-d", "--dirs", help="Individual directories to process", nargs="*")
    group.add_argument("-j", "--trajectories", help="Individual trajectory files to process", nargs="*")

    parser.add_argument("-p", "--plot", action="append", dest="plots", default=[], help="Select plot(s) (available: %s)" % ", ".join(PLOTS))
    parser.add_argument("-m", "--measurement", action="append", dest="measurements", default=[], help="Select measurement(s) (available: %s)" % ", ".join(SimulationSet.MEASUREMENTS))
    parser.add_argument("-c", "--chain", action="append", dest="chains", default=["all"], help="Select chain(s)")
    parser.add_argument("-t", "--temperature", action="append", type=float, dest="temperatures", default=[300.0], help="Select temperature(s)")

    parser.add_argument("-f", "--fit", help="Attempt to fit values to function with the provided exponent (plot vs N only)")

    parser.add_argument("-v", "--verbose", help="Turn on verbose output", action="count")

    args = parser.parse_args()

    if args.verbose > 1:
	logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)

    if not args.dirs and not args.trajectories:
        logging.warn("No directories or trajectories given. Exiting.")
        sys.exit(1)
    elif args.dirs:
        s = SimulationSet.from_dirlist(args)
    else:
        s = SimulationSet.from_trajectorylist(args)
        
    s.all_plots(args)



