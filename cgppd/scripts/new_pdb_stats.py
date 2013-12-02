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


class Residue(object):
    def __init__(self, amino_acid, index, x, y, z):
        self.amino_acid = amino_acid
        self.index = index
        self.position = np.array((x, y, z))

    def __str__(self):
        return "%d: %s %s" % (self.index, self.amino_acid, self.position)


class Chain(object):
    def __init__(self, chain_id, residues=None, length=None, radius=None):
        self.chain_id = chain_id
        self.residues = residues or [] # list
        self.length = length
        self.radius = radius

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
        self.chains = chains or {} # dict

    def add_chain(self, chain):
        self.chains[chain.chain_id] = chain

    def __str__(self):
        return "sample %d at %fK with potential %f" % (self.sample, self.temperature, self.potential)


class Temperature(object):
    def __init__(self, temperature, samples):
        self.temperature = temperature
        self.samples = samples # list

    def __str__(self):
        return "temperature %f with %d samples" % (self.temperature, len(self.samples))


class Simulation(object):
    FILENAME = re.compile('output/(.*)/pdb/sample_(.*)_(.*)K_ *(.*).pdb')
    ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')
    SAMPLE = re.compile('(\d+) at (\d+\.\d+)K with potential (-?\d+\.\d+)\n(.*)', re.MULTILINE | re.DOTALL)
    CHAIN = re.compile('chain (.*) length (\d+\.\d+) radius (\d+\.\d+)')

    def __init__(self, N, temperatures):
        self.N = N
        self.temperatures = temperatures # dict

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
            temperatures = defaultdict(list)

            for s in re.split("sample", summary):
                if not s:
                    continue

                sample, temperature, potential, chains = cls.SAMPLE.search(s).groups()
                proteins = {}

                for c in chains.split("\n"):
                    if not c:
                        continue
                    chain, length, radius = cls.CHAIN.match(c).groups()
                    proteins[chain] = Protein(chain, float(length), float(radius))

                temperatures[float(temperature)].append(Conformation(int(sample), float(temperature), float(potential), proteins))

            # TODO: get chains and name from PDB
            chains = dict(n.split(":") for n in args.names)

            name = SUMMARYFILENAME.match(filename).group(1)

            return cls(N, temperatures)


class SimulationSet(object):
    def __init__(self, simulations):
        self.simulations = simulations # list



