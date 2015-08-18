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

# TODO TODO TODO the sample step and the sample number are not the same thing. Add them both. enumerate files...

class Sample(object):
    def __init__(self, sample_no, length, radius, potential, cluster=None):
        self.sample_no = sample_no
        self.length = length
        self.radius = radius
        self.potential = potential
        self.cluster = cluster # sample_no

    @classmethod
    def from_PDB(cls, pdb_file):
        residues = []
        
        for line in pdb_file:
            if "sample" in line:
                sample_no = int(re.search("(\d+)", line).group(1))
            elif "potential" in line:
                potential = float(re.search("(-?\d+\.\d+)", line).group(1))
            elif "ATOM" in line:
                residues.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

        length, radius = measure(residues)
        
        return cls(sample_no, length, radius, potential)


class Cluster(object):
    LOG_ROW_SEP = re.compile(" *\| *")
    MEMBER_SEP = re.compile(" *")

    def __init__(self, sample_no, length, radius, potential=None, members=None):
        self.sample_no = sample_no
        self.length = length
        self.radius = radius
        self.potential = potential
        self.members = members or [] # sample_nos

    @classmethod
    def from_PDB(cls, pdb_file):
        clusters = {}

        residues = []
        for line in pdb_file:
            if "TITLE" in line:
                sample_no = int(float(re.search("(-?\d+\.\d+)", line).group(1))) + 1
            elif "ATOM" in line:
                residues.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            elif "ENDMDL" in line:
                length, radius = measure(residues)
                clusters[sample_no] = cls(sample_no, length, radius)
                residues = []
                del sample_no
                
        return clusters

    @classmethod
    def match_to_samples(cls, clusters, cluster_log_file, samples):
        for line in cluster_log_file:
            if "|" not in line:
                continue
            elif line.startswith("cl."):
                continue
            else:
                _, _, frame_id, members = cls.LOG_ROW_SEP.split(line)

                if frame_id:
                    if " " in frame_id:
                        frame_id, _ = frame_id.split(" ")

                    cluster = clusters[int(frame_id) + 1]
               
                cluster.members.extend(int(m) + 1 for m in cls.MEMBER_SEP.split(members))

        for cluster in clusters.itervalues():
            cluster.potential = samples[cluster.sample_no].potential
            for sample_no in cluster.members:
                samples[sample_no].cluster = cluster.sample_no


class Simulation(object):
    def __init__(self, samples, clusters):
        self.samples = samples # list
        self.clusters = clusters # dict

    @classmethod
    def from_dir(cls, directory):
        print "Processing directory %s..." % directory
        
        samples = []
        clusters = {}

        summary_filename = os.path.join(directory, "summary.csv")

        if os.path.isfile(summary_filename):
            print "Reading sample summary..."
            with open(summary_filename, "r") as summary_file:
                reader = csv.reader(summary_file)
                for sample_no, length, radius, potential in reader:
                    samples.append(Sample(int(sample_no), float(length), float(radius), float(potential)))
                    
        else:
            print "Writing sample summary..."
                            
            for pdb_filename in glob.iglob(os.path.join(directory, "pdb", "sample*_300.0K_*.pdb")):
                with open(pdb_filename, "r") as pdbfile:
                    samples.append(Sample.from_PDB(pdbfile))

            samples.sort(key=attrgetter("sample_no"))

            with open(summary_filename, "w") as summary_file:
                writer = csv.writer(summary_file)
                for sample in samples:
                    writer.writerow([sample.sample_no, sample.length, sample.radius, sample.potential])

        # clusters
        
        cluster_summary_filename = os.path.join(directory, "cluster_summary.csv")
        cluster_log_filename = os.path.join(directory, "cluster.log")
        cluster_filename = os.path.join(directory, "clusters.pdb")

        if os.path.isfile(cluster_summary_filename):
            print "Reading cluster summary..."
            with open(cluster_summary_filename, "r") as cluster_summary_file:
                reader = csv.reader(cluster_summary_file)
                for sample_no, length, radius, potential, members_str in reader:
                    members = [int(m) for m in members_str.split()]
                    
                    cluster = Cluster(int(sample_no), float(length), float(radius), float(potential), members)
                    clusters.append(cluster)
                    
                    for m in members:
                        samples[m].cluster = cluster.sample_no
                        
        elif os.path.isfile(cluster_filename) and os.path.isfile(cluster_log_filename):
            print "Writing cluster summary..."
            with open(cluster_filename, "r") as cluster_file:
                clusters = Cluster.from_PDB(cluster_file) # without potential or members

            with open(cluster_log_filename, "r") as cluster_log_file:
                Cluster.match_to_samples(clusters, cluster_log_file, samples)

            with open(cluster_summary_filename, "w") as cluster_summary_file:
                writer = csv.writer(cluster_summary_file)

                for sample_no, cluster in sorted(clusters.iteritems()):
                    writer.writerow([sample_no, cluster.length, cluster.radius, cluster.potential, " ".join(str(m) for m in cluster.members)])
                
        else:
            print "No cluster information found. Please run g_cluster."

        return cls(samples, clusters)

def plot_length(simulations):
    print "Plotting height"
    for simulation in simulations:
        plt.hist([s.length for s in simulation.samples])
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("dirs", help="Individual directories to process", nargs="+")

    args = parser.parse_args()

    simulations = [Simulation.from_dir(d) for d in args.dirs]
    plot_length(simulations)

