#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
from collections import defaultdict

COMMENT = re.compile(r"^\s*#")

def pair_potentials_from_file(filename):
    with open(filename, 'r') as pp_file:
        lines = [l for l in pp_file if not COMMENT.search(l)]

    amino_acids = [a.upper() for a in lines[0].strip().split()]

    pair_potentials = defaultdict(dict)

    for l in lines[1:]:
        chunks = l.strip().split()
        acid_i = chunks[0].upper()
        other_acids = amino_acids[amino_acids.index(acid_i):]
        potentials = [float(p) for p in chunks[1:]]

        for acid_j, potential_value in zip(other_acids, potentials):
            pair_potentials[acid_i][acid_j] = potential_value
            pair_potentials[acid_j][acid_i] = potential_value
            print "%s-%s: %g" % (acid_i, acid_j, potential_value)

    return pair_potentials


#class Potential(object):

    #def __init__(self, residues, pair_potentials):
        #self.residues = residues
        #self.pair_potentials = pair_potentials
        #self.components = {}

    #@property
    #def LJ(self):
        #if "LJ" not in self.components:
            #pass # TODO

        #return self.components["LJ"]

    #@property
    #def DH(self):
        #pass # TODO

    #@property
    #def bond(self):
        #pass # TODO

    #@property
    #def angle(self):
        #pass # TODO

    #@property
    #def torsion(self):
        #pass # TODO

    #@property
    #def total(self):
        #pass # TODO

if __name__ == "__main__":
     potentials = pair_potentials_from_file(sys.argv[1])
