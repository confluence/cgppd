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
import csv

import logging

class BananaError(Exception):
    pass


class Residue(object):
    def __init__(self, amino_acid, index, x, y, z):
        self.amino_acid = amino_acid
        self.index = index
        self.position = np.array((x, y, z))


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


class Sample(object):
    def __init__(self, sample_id, temperature, potential, chains=None, cluster=None):
        self.sample_id = sample_id
        self.temperature = temperature
        self.potential = potential
        self.chains = chains or []
        self.cluster = cluster

    def add_chain(self, chain):
        self.chains.append(chain)


class Temperature(object):
    def __init__(self, temperature, samples=None, cluster_set=None):
        self.temperature = temperature
        self.samples = samples or []
        self.cluster_set = cluster_set or []

    def add_sample(self, sample):
        self.samples.append(sample)

    def add_clusters(self, cluster_set):
        self.cluster_set = cluster_set
        for cluster in cluster_set.clusters.itervalues():
            cluster.match_to_samples(self.samples)


class Cluster(object):
    def __init__(self, frame_id, chains=None, members=None):
        self.frame_id = frame_id
        self.chains = chains or []
        self.members = members or []

    def add_chain(self, chain):
        self.chains.append(chain)

    def match_to_samples(self, samples):
        self.samples = [samples[i + 1] for i in self.members]


class ClusterSet(object): 
    FRAME = re.compile('TITLE *frame t= (-?\d+).000')
    LOG_ROW_SEP = re.compile(" *\| *")
    MEMBER_SEP = re.compile(" *")

    def __init__(self, clusters):
        self.clusters = clusters # dict with frame as key
    
    @classmethod
    def write_summary(cls, trajectoryfilename, logfilename, clustersummaryfilename, args):
        with open(trajectoryfilename, 'r') as trajectoryfile:
            clusters = {}
            logging.info("Parsing trajectory file '%s' ..." % (trajectoryfilename))
            
            for line in trajectoryfile:
                if "frame" in line: # clustered trajectory
                    cluster = Cluster(int(cls.FRAME.match(line).group(1)))
                    chain = None
                elif line.startswith("ATOM"):
                    amino_acid, chain_id, index, x, y, z = Simulation.ATOM.match(line).groups()

                    if not chain:
                        chain = Chain(chain_id)
                    chain.add_residue(Residue(amino_acid, int(index), float(x), float(y), float(z)))

                elif line.startswith("TER"): # end of chain
                    cluster.add_chain(chain)
                    chain = None
                
                elif line.startswith("ENDMDL"): # end of cluster
                    # quick and dirty hack to calculate stats for multi-chain molecule
                    # length relies on sensible ordering of chains in the PDB file
                    molecule = Chain("all")
                    
                    for chain in cluster.chains:
                        chain.measure()
                        molecule.residues.extend(chain.residues)

                    molecule.measure()
                    cluster.add_chain(molecule)
                    clusters[cluster.frame_id] = cluster

        with open(logfilename, 'r') as logfile:
            cluster = None;
            
            for line in logfile:
                if "|" not in line:
                    continue
                elif line.startswith("cl."):
                    continue
                else:
                    _, _, frame_id, members = cls.LOG_ROW_SEP.split(line)

                    if frame_id:
                        if " " in frame_id:
                            frame_id, _ = frame_id.split(" ")
                        cluster = clusters[int(frame_id)]
                    
                    cluster.members.extend(int(m) for m in cls.MEMBER_SEP.split(members))

        with open(clustersummaryfilename, 'w') as summaryfile:
            writer = csv.writer(summaryfile)
            for frame_id, cluster in clusters.iteritems():
                row = [frame_id, " ".join(str(m) for m in cluster.members)]
                
                # quick and dirty hack to calculate stats for multi-chain molecule
                # length relies on sensible ordering of chains in the PDB file
                molecule = Chain("all")
                
                for chain in cluster.chains:
                    chain.measure()
                    row.extend([chain.chain_id, chain.length, chain.radius])

                    molecule.residues.extend(chain.residues)

                molecule.measure()
                row.extend([molecule.chain_id, molecule.length, molecule.radius])

                writer.writerow(row)

    @classmethod
    def from_summary(cls, filename, args):
        clusters = {}
        
        with open(filename, 'r') as summaryfile:
            reader = csv.reader(summaryfile)
            for row in reader:
                frame_id, members = row[:2]
                cluster = Cluster(int(frame_id), members=[int(m) for m in members.split()])

                for i in range(2, len(row), 3):
                    chain_id, length, radius = row[i:i+3]
                    cluster.add_chain(Chain(chain_id, float(length), float(radius)))

                clusters[cluster.frame_id] = cluster

        return cls(clusters)


class Simulation(object):
    FILENAME = re.compile('.*/pdb/sample_(.*)_(.*)K_ *(.*).pdb')
    NAME = re.compile('(\d*)$')

    ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')

    def __init__(self, temperatures, name, N=None):
        self.temperatures = temperatures # actually needs to be a dict
        self.name = name
        self.N = N

    @classmethod
    def write_summary(cls, filenames, summaryfilename, args):
        num_files = len(filenames)
        with open(summaryfilename, 'w') as summaryfile:
            writer = csv.writer(summaryfile)
            
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

                    row = [sample.sample_id, sample.temperature, sample.potential]
                    
                    ## quick and dirty hack to calculate stats for multi-chain molecule
                    ## length relies on sensible ordering of chains in the PDB file
                    molecule = Chain("all")
                    
                    for chain in sample.chains:
                        chain.measure()
                        row.extend([chain.chain_id, chain.length, chain.radius])

                        molecule.residues.extend(chain.residues)

                    molecule.measure()
                    row.extend([molecule.chain_id, molecule.length, molecule.radius])
                    writer.writerow(row)

    @classmethod
    def from_summary(cls, filename, args):
        directory, _ = os.path.split(filename)
        
        rest, last_dir = os.path.split(directory)
        name = "_".join(last_dir.split("_")[:-1])

        N = name[3:] if name.startswith("ala") else None

        #N = cls.NAME.search(name).group(1) or None

        temperatures = {}

        with open(filename) as summaryfile:
            reader = csv.reader(summaryfile)

            for row in reader:
                sample_id, temperature, potential = row[:3]
                sample = Sample(int(sample_id), float(temperature), float(potential))

                for i in range(3, len(row), 3):
                    chain_id, length, radius = row[i:i+3]
                    sample.add_chain(Chain(chain_id, float(length), float(radius)))

                if sample.temperature not in temperatures:
                    temperatures[sample.temperature] = Temperature(sample.temperature)

                temperatures[sample.temperature].add_sample(sample)

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
    def from_dirlist(cls, args):
        simulations = []

        for d in args.dirs:
            logging.info("Processing directory '%s'..." % d)

            summaryfilename = "%s/summary.csv" % d
            sim_dir = "%s/pdb" % d

            if not os.path.isfile(summaryfilename):
                logging.info("Writing new summary file...")
                if os.path.exists(sim_dir):
                    files = glob.glob("%s/*" % sim_dir)
                    Simulation.write_summary(files, summaryfilename, args)

            logging.info("Loading summary file...")
            s = Simulation.from_summary(summaryfilename, args)

            if 300.0 in s.temperatures:
                clusterfilename = os.path.join(d, "clusters.pdb")
                clusterlogname = os.path.join(d, "cluster.log")

                if os.path.isfile(clusterfilename) and os.path.isfile(clusterlogname):
                    logging.info("Cluster data detected.")
                    clustersummaryfilename = "%s/cluster_summary.csv" % d

                    if not os.path.isfile(clustersummaryfilename):
                        logging.info("Writing new cluster summary file...")                    
                        ClusterSet.write_summary(clusterfilename, clusterlogname, clustersummaryfilename, args)

                    logging.info("Loading cluster summary file...")
                    s.temperatures[300.0].add_clusters(ClusterSet.from_summary(clustersummaryfilename, args))

            simulations.append(s)            

        return cls(simulations)

    def _display_plot(self, plt, fig, title, xlabel, ylabel):
        plt.tight_layout()
        fig.suptitle(title, fontsize=14)
        fig.text(0.5, 0.02, xlabel, ha='center', va='center')
        fig.text(0.02, 0.5, ylabel, ha='center', va='center', rotation='vertical')
        plt.show()
        plt.close()

    def plot_cluster_members(self, dataset, args):
        if 300.0 not in dataset.temperatures:
            raise BananaError("Could not find 300K in available temperatures.")

        pass # TODO; what should we actually plot? cluster vs % of population?


    # TODO plot all simulations on the same axes
            

    def plot_histogram(self, dataset, args):
        print dataset.simulations
        for measurement in dataset.measurements:

            for sim in dataset.simulations:
                print sim
                for chain_id in dataset.chain_ids:
                    fig = plt.figure()

                    for i, temperature in enumerate(dataset.temperatures):
                        ax = fig.add_subplot(dataset.dimy, dataset.dimx, i + 1)

                        chains = dataset.data[(sim, chain_id, temperature)]

                        values = np.array([getattr(c, measurement) for c in chains])
                        ax.hist(values, bins=50)

                        ax.set_title("%sK" % temperature)

                        labels = ax.get_xticklabels()
                        for label in labels:
                            label.set_rotation(30)
                            label.set_fontsize(10)

                    title = "%s %s %s" % (sim, chain_id, measurement)
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

                    title = "%s %s %s" % (sim, chain_id, measurement)
                    xlabel = "Sample time"
                    ylabel = u"%s (Å)" % measurement

                    self._display_plot(plt, fig, title, xlabel, ylabel)

    def plot_vs_N(self, dataset, args):
        try:
            xvalues = np.array([int(N) for N in dataset.simulations])
        except TypeError:
            raise BananaError("This dataset does not have a range of residue sizes.")

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
                continue

            plot_method(dataset, args)


PLOTS = tuple(n[5:] for n in SimulationSet.__dict__ if n.startswith("plot_"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")

    parser.add_argument("dirs", help="Individual directories to process", nargs="+")

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

    try:
        s = SimulationSet.from_dirlist(args)
        s.all_plots(args)
    except BananaError, e:
        logging.error(e)
        sys.exit(1)



