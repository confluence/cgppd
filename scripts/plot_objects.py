#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import csv
import sys
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
    def __init__(self, sample_no, sample_step, length, radius, potential, cluster=None):
        self.sample_no = sample_no
        self.sample_step = sample_step
        self.length = length
        self.radius = radius
        self.potential = potential
        self.cluster = cluster # cluster object

    @classmethod
    def from_PDB(cls, pdb_file):
        residues = []
        
        for line in pdb_file:
            if "sample" in line:
                sample_step = int(re.search("(\d+)", line).group(1))
            elif "potential" in line:
                potential = float(re.search("(-?\d+\.\d+)", line).group(1))
            elif "ATOM" in line:
                # we can deal with overflows here which are technically illegal PDB syntax, because we know the precision is fixed
                residues.append(tuple(float(n) for n in re.findall('(-?\d+\.\d{3})', line[30:-12])))

        length, radius = measure(residues)
        
        return cls(None, sample_step, length, radius, potential)


class Cluster(object):
    LOG_ROW_SEP = re.compile(" *\| *")
    MEMBER_SEP = re.compile(" *")

    def __init__(self, sample_no, length, radius, potential=None, member_nos=None, members=None):
        self.sample_no = sample_no
        self.length = length
        self.radius = radius
        self.potential = potential
        self.member_nos = member_nos or [] # sample_nos
        self.members = members or [] # sample objects

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
               
                cluster.member_nos.extend(int(m) + 1 for m in cls.MEMBER_SEP.split(members))

        for cluster in clusters.itervalues():
            cluster.potential = samples[cluster.sample_no].potential
            for sample_no in cluster.member_nos:
                samples[sample_no].cluster = cluster
                cluster.members.append(samples[sample_no])


class Simulation(object):
    def __init__(self, samples, clusters):
        self.samples = samples # list
        self.clusters = clusters # dict

    @classmethod
    def from_dir(cls, directory, cutoff=None):
        print "Processing directory %s..." % directory
        
        samples = []
        clusters = {}

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
            print "Detecting temperature nearest 300K..."

            # we are being as agnostic as possible about the number of replicas, steps before first sample, etc..
            temps_seen = set()
            closest_temp = None

            # sorting is horrible, but the only way
            for pdb_filename in sorted(glob.iglob(os.path.join(directory, "pdb", "sample*.pdb"))):
                temp = float(re.search("_(\d{3}\.\d)K_*", pdb_filename).group(1))
                
                if temp in temps_seen:
                    break
                    
                if closest_temp is None or abs(300 - temp) < abs(300 - closest_temp):
                    closest_temp = temp

                temps_seen.add(temp)

            print "Using %.1fK." % closest_temp
            
            print "Writing sample summary..."
                            
            for pdb_filename in glob.iglob(os.path.join(directory, "pdb", "sample*_%.1fK_*.pdb" % closest_temp)):
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

        # clusters

        if cutoff:
            print "Using cluster with cutoff %s" % cutoff
            cutoff_dirname = "clusters_%s" % cutoff
        
            cluster_summary_filename = os.path.join(directory, cutoff_dirname, "cluster_summary.csv")
            cluster_header = ("sample_no", "length", "radius", "potential", "members")
            cluster_log_filename = os.path.join(directory, cutoff_dirname, "cluster.log")
            cluster_filename = os.path.join(directory, cutoff_dirname, "clusters.pdb")

            if os.path.isfile(cluster_summary_filename):
                print "Reading cluster summary..."
                with open(cluster_summary_filename, "r") as cluster_summary_file:
                    reader = csv.reader(cluster_summary_file)
                    reader.next() # skip header
                    for sample_no, length, radius, potential, members_str in reader:
                        member_nos = [int(m) for m in members_str.split()]
                        
                        cluster = Cluster(int(sample_no), float(length), float(radius), float(potential), member_nos)
                        clusters[int(sample_no)] = cluster
                        
                        for m in member_nos:
                            samples[m].cluster = cluster
                            cluster.members.append(samples[m])
                            
            elif os.path.isfile(cluster_filename) and os.path.isfile(cluster_log_filename):
                print "Writing cluster summary..."
                with open(cluster_filename, "r") as cluster_file:
                    clusters = Cluster.from_PDB(cluster_file) # without potential or members

                with open(cluster_log_filename, "r") as cluster_log_file:
                    Cluster.match_to_samples(clusters, cluster_log_file, samples)

                with open(cluster_summary_filename, "w") as cluster_summary_file:
                    writer = csv.writer(cluster_summary_file)
                    writer.writerow(cluster_header)

                    for sample_no, cluster in sorted(clusters.iteritems()):
                        writer.writerow([sample_no, cluster.length, cluster.radius, cluster.potential, " ".join(str(m) for m in cluster.members)])
                
            else:
                print "No cluster information found."

        return cls(samples, clusters)
