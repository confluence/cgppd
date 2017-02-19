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
import h5py
from collections import defaultdict

def measure(residues):
    length = np.linalg.norm(np.array(residues[-1]) - np.array(residues[0]))

    mass = len(residues)
    centre_of_mass = np.array(residues).mean(axis=0)
    # Avoid having to square unnecessary square roots by calculating squared distances with a dot product
    diff_vectors = [(np.array(r) - centre_of_mass) for r in residues]
    radius = np.sqrt(sum(np.dot(v.T, v) for v in diff_vectors) / mass)

    return length, radius

class Sample(object):
    XYZ = re.compile('(-?\d+\.\d{3})')
    
    def __init__(self, sample_no, sample_step, length, radius, potential):
        self.sample_no = sample_no
        self.sample_step = sample_step
        self.length = length
        self.radius = radius
        self.potential = potential

    @classmethod
    def from_PDB(cls, pdb_file):
        samples = []
        
        residues = []
        sample_step=None
        potential = None
        
        for line in pdb_file:
            if "sample" in line:
                sample_step = int(re.search("(\d+)", line).group(1))
            elif "potential" in line:
                potential = float(re.search("(-?\d+\.\d+)", line).group(1))
            elif "ATOM" in line:
                # we can deal with overflows here which are technically illegal PDB syntax, because we know the precision is fixed
                residues.append(tuple(float(n) for n in cls.XYZ.findall(line[30:-12])))
            elif "END" in line:
                length, radius = measure(residues)
                samples.append(cls(None, sample_step, length, radius, potential))
                residues = []
                sample_step=None
                potential = None

        return samples
        

class Cluster(object):
    def __init__(self, sample_nos):
        self.sample_nos = sorted(sample_nos)
        self.samples = None
        
    def match_to_samples(self, samples):
        self.samples = [samples[i] for i in self.sample_nos]
        
    @classmethod
    def from_string(cls, sample_nos_string):
        return cls([int(i) for i in sample_nos_string.split()])
        
        
class SampleContacts(object):
    XYZ = re.compile('(-?\d+\.\d{3})')
    
    def __init__(self, contact_filename):
        self.contact_filename = contact_filename
        
    def average_contacts(self, cutoff):
        f = h5py.File(self.contact_filename, "r")
        
        contacts = []
        
        for sample in sorted(f.keys()):
            distances = f["sample"]
            within_cutoff = (-1 < distances < cutoff).sum(axis=1)
            contacts.append(within_cutoff)
            
        contacts = np.array(contacts)
        averages = contacts.mean(axis=1)
        
        f.close()
        return averages
                
    def write_records_from_pdb(self, pdbfile, sample_no=0):
        samples = []
        residues = []
        
        for line in pdb_file:
            if "ATOM" in line:
                chain = line[21]
                # we can deal with overflows here which are technically illegal PDB syntax, because we know the precision is fixed
                pos = np.array([float(n) for n in self.XYZ.findall(line[30:-12])])
                residues.append((chain, pos))
            elif "END" in line:
                distances = []
                for (chain1, pos1) in residues:
                    row = []
                    for (chain2, pos2) in residues:
                        if chain1 == chain2:
                            row.append(-1)
                        else:
                            row.append(np.linalg.norm(pos1 - pos2))
                    distances.append(row)
                
                samples.append(distances)
                residues = []

        f = h5py.File(self.contact_filename, "w")
        
        for i, distances in enumerate(samples, sample_no):
            f.create_dataset("sample_%08d" % i, data=np.array(distances))
                
        f.close()


class Simulation(object):    
    def __init__(self, samples, cluster_sets, contacts):
        self.samples = samples # list
        self.cluster_sets = cluster_sets # dict
        self.contacts = contacts

    @classmethod
    def from_dir(cls, directory, args, calculate_contacts=True):
        print "Processing directory %s..." % directory
        
        samples = []

        summary_filename = os.path.join(directory, "summary.csv")
        header = ("sample_no", "sample_step", "length", "radius", "potential")
        
        trajectory_filename = os.path.join(directory, "trajectory.pdb")
        
        write = False
        
        def get_closest_temp():
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
            return closest_temp

        if os.path.isfile(summary_filename):
            print "Reading sample summary..."
            with open(summary_filename, "r") as summary_file:
                reader = csv.reader(summary_file)
                reader.next() # skip header
                for sample_no, sample_step, length, radius, potential in reader:
                    samples.append(Sample(int(sample_no), int(sample_step or "0"), float(length), float(radius), float(potential or "0")))
                    
        elif os.path.isfile(trajectory_filename):
            print "Using trajectory file..."
            
            with open(trajectory_filename, "r") as pdbfile:
                samples.extend(Sample.from_PDB(pdbfile))
                
            write = True
        else:
            closest_temp = get_closest_temp()
                            
            for pdb_filename in glob.iglob(os.path.join(directory, "pdb", "sample*_%.1fK_*.pdb" % closest_temp)):
                with open(pdb_filename, "r") as pdbfile:
                    samples.extend(Sample.from_PDB(pdbfile))
                    
            samples.sort(key=attrgetter("sample_step"))
                    
            write = True
            
        if write:
            print "Writing sample summary..."
            
            for i, sample in enumerate(samples):
                sample.sample_no = i

            with open(summary_filename, "w") as summary_file:
                writer = csv.writer(summary_file)
                writer.writerow(header)
                for sample in samples:
                    writer.writerow([sample.sample_no, sample.sample_step, sample.length, sample.radius, sample.potential])
                    
        # new cluster stuff
                
        cluster_sets = {}
        
        cluster_filenames = glob.glob(os.path.join(directory, "clusters*.txt"))
        
        for cluster_filename in cluster_filenames:
            cluster_description, _ = os.path.splitext(os.path.basename(cluster_filename))
            clusters = []

            with open(cluster_filename, "r") as cluster_file:
                for line in cluster_file:
                    cluster = Cluster.from_string(line)
                    cluster.match_to_samples(samples)
                    clusters.append(cluster)
            
            cluster_sets[cluster_description] = clusters
            
        if not cluster_filenames:
            print "No cluster information found."
            
        # Contact averages
        
        if calculate_contacts:
            contact_filename = os.path.join(directory, "contacts.hdf5")
            contacts = Contacts(contact_filename)
                
            if os.path.isfile(contact_filename):
                print "Reading contacts from file..." # That's it; the file will be opened by the object
            else:
                if os.path.isfile(trajectory_filename):
                    print "Using trajectory file for contacts..."
                
                    with open(trajectory_filename, "r") as pdbfile:
                        contacts.write_records_from_pdb(pdbfile)
                        
                else:
                    print "Using sample files for contacts..."
                    closest_temp = get_closest_temp()
                                    
                    for sample_no, pdb_filename in enumerate(glob.iglob(os.path.join(directory, "pdb", "sample*_%.1fK_*.pdb" % closest_temp))):
                        with open(pdb_filename, "r") as pdbfile:
                            contacts.write_records_from_pdb(pdbfile, sample_no)
        else:
            contacts = None

        return cls(samples, cluster_sets, directory, contacts)
    


class PolyalanineSimulationSequence(object):
    NAME = re.compile("ala(\d+)_.*")
    
    def __init__(self, sims):
        self.sims = sims
        
    @classmethod
    def from_dirs(cls, args):
        sims = []
        
        for d in args.dirs:
            dirname = os.path.basename(d)
            name_match = cls.NAME.match(dirname)
            if name_match is None:
                sys.exit("'%s' does not look like a polyalanine simulation." % dirname)
            n = int(name_match.group(1))
            
            sims.append((n, Simulation.from_dir(d, args, calculate_contacts=False)))
            
        return cls(sorted(sims))
        
        
class DiubiquitinSimulationGroup(object):
    NAME = re.compile("diubiquitin_(lys|met)_(\d+)_.*")
    
    def __init__(self, sims):
        self.sims = sims
        self.extra_properties = set()
        
    @classmethod
    def from_dirs(cls, args):
        sims = []
        
        for d in args.dirs:
            dirname = os.path.basename(d)
            name_match = cls.NAME.match(dirname)
            if name_match is None:
                sys.exit("'%s' does not look like a diubiquitin simulation." % dirname)
            res, index = name_match.groups()
            
            sims.append(("%s-%s" % (res.upper(), index), Simulation.from_dir(d, args)))
            
        return cls(sims)
