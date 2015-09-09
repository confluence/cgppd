#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import csv
import os
import glob
import argparse
from operator import attrgetter
import numpy as np
import matplotlib.pyplot as plt

def measure(residues):
    length = np.linalg.norm(np.array(residues[-1]) - np.array(residues[0]))

    mass = len(residues)
    centre_of_mass = np.array(residues).mean(axis=0)
    # Avoid having to square unnecessary square roots by calculating squared distances with a dot product
    diff_vectors = [(np.array(r) - centre_of_mass) for r in residues]
    radius = np.sqrt(sum(np.dot(v.T, v) for v in diff_vectors) / mass)

    return length, radius


class Sample(object):
    def __init__(self, sample_no, sample_step, length, radius, potential):
        self.sample_no = sample_no
        self.sample_step = sample_step
        self.length = length
        self.radius = radius
        self.potential = potential

    @classmethod
    def from_PDB(cls, pdb_file):
        residues = []
        
        for line in pdb_file:
            if "sample" in line:
                sample_step = int(re.search("(\d+)", line).group(1))
            elif "potential" in line:
                potential = float(re.search("(-?\d+\.\d+)", line).group(1))
            elif "ATOM" in line:
                residues.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

        length, radius = measure(residues)
        
        return cls(None, sample_step, length, radius, potential)

class Simulation(object):
    def __init__(self, samples):
        self.samples = samples # list

    @classmethod
    def from_dir(cls, directory):
        print "Processing directory %s..." % directory
        
        samples = []

        summary_filename = os.path.join(directory, "summary.csv")
        header = ("sample_no", "sample_step", "length", "radius", "potential")

        if os.path.isfile(summary_filename):
            print "Reading sample summary..."
            with open(summary_filename, "r") as summary_file:
                reader = csv.reader(summary_file)
                reader.next() # skip header
                for sample_no, sample_step, length, radius, potential in reader:
                    samples.append(Sample(int(sample_no), int(sample_step), float(length), float(radius), float(potential)))
                    
        else:
            print "Writing sample summary..."
                            
            for pdb_filename in glob.iglob(os.path.join(directory, "pdb", "sample*_300.0K_*.pdb")):
                with open(pdb_filename, "r") as pdbfile:
                    samples.append(Sample.from_PDB(pdbfile))

            samples.sort(key=attrgetter("sample_step"))
            
            for i, sample in enumerate(samples):
                sample.sample_no = i

            with open(summary_filename, "w") as summary_file:
                writer = csv.writer(summary_file)
                writer.writerow(header)
                for sample in samples:
                    writer.writerow([sample.sample_no, sample.sample_step, sample.length, sample.radius, sample.potential])

        return cls(samples)


class PolyalanineSimulationSequence(object):
    NAME = re.compile("ala(\d+)_.*")
    
    def __init__(self, sims):
        self.sims = sims
        
    @classmethod
    def from_dirs(cls, dirs):
        sims = []
        
        for d in dirs:
            dirname = os.path.basename(d)
            name_match = NAME.match(dirname)
            if name_match is None:
                sys.exit("'%s' does not look like a polyalanine simulation." % dirname)
            n = int(name_match.group(1))
            
            sims.append((n, Simulation.from_dir(d)))
            
        return cls(sorted(sims))

# TODO: use a simulation set; commandline param for type of LJ potential
# TODO: refactor some of this; maybe factor out classes into a library

def plot_mean_radius(simulation_set, lj):
    xvalues, sims = zip(*simulation_set.sims)
    values = [np.sqrt(np.mean(s.radius**2)) for sim in sims for s in sim.samples]
    
    plt.plot(xvalues, values, 'bo')
    
    if lj == "off":
        power = 1.0/2.0
    elif lj == "repulsive":
        power = 2.0/3.0

    fitfunc = lambda p, N: p[0] * np.array(N) ** power # Target function
    errfunc = lambda p, N, y: fitfunc(p, N) - y # Distance to the target function
    p0 = [1.0] # Initial guess for the parameters
    p1, success = optimize.leastsq(errfunc, p0[:], args=(xvalues, values))

    #logging.debug("A = %g" % p1[0])

    plt.plot(xvalues, [fitfunc(p1, x) for x in xvalues], 'b-')

    plt.set_title("%sK" % temperature)


    plt.set_xscale('log', basex=2) # TODO: investigate using the loglog function instead; maybe add an option for it
    plt.set_yscale('log')
        
    plt.show()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("lj", help="Type of Lennard-Jones interaction", choices=("off", "repulsive"))
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")

    args = parser.parse_args()

    simulation_set = PolyalanineSimulationSequence.from_dirlist(args.dirs)
    plot_mean_radius(simulation_set, args.lj)
