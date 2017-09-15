#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import csv
import sys
import os
import glob
from operator import attrgetter
import numpy as np
from collections import defaultdict
import h5py

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
        
        
class Contacts(object):
    XYZ = re.compile('(-?\d+\.\d{3})')
    
    def __init__(self, contact_filename, num_samples):
        self.contact_filename = contact_filename
        self.num_samples = num_samples
        
    def average_contacts(self, cutoff):
        within_cutoff = defaultdict(int)
        
        f = h5py.File(self.contact_filename, "r")
        
        
        lhs_residues = [k for k in f.keys() if k != "chains"]
        
        for i in lhs_residues:
            for j in f[i].keys():
                samples_within_cutoff = sum(np.array(f["%s/%s" % (i, j)]["distances"]) < cutoff)
                within_cutoff[int(i)] += samples_within_cutoff
                within_cutoff[int(j)] += samples_within_cutoff
        
        chains = f["chains"].keys()
        
        averages_by_chain = []
        
        for chain in chains:
            chain_residues = f["chains"][chain]["residues"]
            averages = [within_cutoff[i] / float(self.num_samples) for i in chain_residues]
            averages_by_chain.append((chain, averages))
        
        f.close()
        return averages_by_chain

    @staticmethod
    def cached_contacts_from_file(filename):
        print "Reading cached contact averages..."
        
        averages_by_chain = {}
        
        with open(filename, "r") as f:
            reader = csv.reader(f)
            reader.next() # skip header
            for chain, averages in reader:
                averages = [float(a) for a in averages.split()]
                averages_by_chain[chain] = averages

        return averages_by_chain.items()

    @staticmethod
    def write_cached_contacts_to_file(filename, averages_by_chain):
        print "Writing cached contact averages..."
                
        with open(filename, "w") as f:
            writer = csv.writer(f)

            writer.writerow(["chain", "resid", "within_cutoff"])

            for chain, averages in averages_by_chain:
                averages = " ".join([str(a) for a in averages])
                writer.writerow([chain, averages])
                
    def write_records_from_pdb(self, pdb_file, sample_no=0):
        f = h5py.File(self.contact_filename, "w")
        
        residues = []
        
        distance_chunks = defaultdict(list)
        BUFFER_SIZE = 10000 # we actually have enough memory and it's faster
        offset = 0
        
        for line in pdb_file:
            if "ATOM" in line:
                chain = line[21]
                # we can deal with overflows here which are technically illegal PDB syntax, because we know the precision is fixed
                pos = np.array([float(n) for n in self.XYZ.findall(line[30:-12])])
                residues.append((chain, pos))
            elif "END" in line:
                for i, (chain1, pos1) in enumerate(residues):
                    for j, (chain2, pos2) in enumerate(residues):
                        if i < j and chain1 != chain2:
                            distance = np.linalg.norm(pos1 - pos2)
                            distance_chunks[(i,j)].append(distance)
                            
                if not (sample_no + 1) % BUFFER_SIZE or sample_no == (self.num_samples - 1):
                    for (i, j), distances in distance_chunks.items():
                        pair = f.require_group("%d/%d" % (i, j))
                        sample_distances = pair.require_dataset("distances", (self.num_samples,), distance.dtype)
                        sample_distances[offset:(sample_no + 1)] = np.array(distances)
                    print "Wrote distances for samples %d to %d..." % (offset, sample_no)
                    distance_chunks = defaultdict(list)
                    offset += BUFFER_SIZE
                            
                    if "chains" not in f:
                        print "Writing chain information..."
                        chains = set(c for c, p in residues)
                        for chain in chains:
                            chain_residues = np.array([i for i, (c, p) in enumerate(residues) if c == chain])
                            chain_grp = f.create_group("chains/%s" % chain)
                            chain_grp.create_dataset("residues", data=chain_residues)
            
                residues = []
                sample_no += 1
                    
        f.close()


class Simulation(object):    
    def __init__(self, samples, cluster_sets, contacts=None, cached_contacts=None):
        self.samples = samples # list
        self.cluster_sets = cluster_sets # dict
        
        self.contacts = contacts
        self.cached_contacts = cached_contacts


    @classmethod
    def from_dir(cls, directory, args):
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

        # Average contacts
        
        contacts = None
        cached_contacts = None

        if hasattr(args, "contact_cutoff"):
            cutoff = args.contact_cutoff

            # Find cached values for this cutoff
            cached_contact_filename = os.path.join(directory, "contacts_%g.csv" % cutoff)

            if os.path.isfile(cached_contact_filename):
                cached_contacts = Contacts.cached_contacts_from_file(cached_contact_filename)

            else:
                # Otherwise load the ginormous full contacts file and calculate this
                contact_filename = os.path.join(directory, "contacts.hdf5")
                contacts = Contacts(contact_filename, len(samples))
                    
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
                                        
                        for i, pdb_filename in enumerate(glob.iglob(os.path.join(directory, "pdb", "sample*_%.1fK_*.pdb" % closest_temp))):
                            with open(pdb_filename, "r") as pdbfile:
                                contacts.write_records_from_pdb(pdbfile, i)

                # cache this cutoff in a file

                cached_contacts = contacts.average_contacts(cutoff)
                Contacts.write_cached_contacts_to_file(cached_contact_filename, cached_contacts)
            
        return cls(samples, cluster_sets, contacts, cached_contacts)
    


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
            
            sims.append((n, Simulation.from_dir(d, args)))
            
        return cls(sorted(sims))
        
        
class DiubiquitinSimulationGroup(object):
    NAME = re.compile("(di|tetra|octa)ubiquitin_(lys|met)_(\d+)_.*")
    
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
            num, res, index = name_match.groups()
            
            sims.append(("%s%s-linked %subq" % (res.title(), index, num), Simulation.from_dir(d, args)))
            
        return cls(sims)
