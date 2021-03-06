#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
from collections import defaultdict
import numpy as np
import argparse

from segment_mockup import Graph

E0 = -2.27
LAMBDA = 0.159
XI = 10.0
E_CHARGE = 1.602176487e-19
EPS = 1e-38

T_ROOM = 298.0
RGAS = 8.3144621
KCAL = 4184.0
N_A = 6.02214129e23
EPS_0 = 8.85418782e-12
ANG = 1e-10
D_WATER = 80.0

RT_TO_KCALMOL = (T_ROOM * RGAS / KCAL) # conversion from RT units to kcal/mol
DH_CONVERSION_FACTOR = (E_CHARGE**2 * N_A / (4.0 * np.pi * EPS_0 * KCAL * ANG * D_WATER)) # as CCELEC in CHARMM

K_SPRING = 378.0
R0 = 3.81

GAMMA_ANGLE =  0.1
EPSILON_ALPHA = 4.3
THETA_ALPHA = 1.6
THETA_BETA = 2.27
K_ALPHA = 106.4
K_BETA = 26.3

COMMENT = re.compile(r"^\s*#")
ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z]?) *\d+ *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) .*')
AMINOACID = re.compile("""#Amino acid class initialisation data
#name
#abbriviation
#short symbol
#van der Waal radius
#residue charge
.*
(.*)
.*
(.*)
(.*)
""")
TORSION = re.compile(r"([A-Z]{3}) +([A-Z]{3}) +(.*) +(\d) +(.*)")

class Residue(object):
    def __init__(self, amino_acid, x, y, z, chain):
        self.amino_acid = amino_acid
        self.position = np.array((x, y, z))
        self.chain = chain
        self.rigid_domain = None
        self.segment_bond = None


class Protein(object):
    def __init__(self):
        self.residues = []
        self.segments = [] # flexible segments as residue indices
        self.all_flexible = False
        self.graph = Graph() # graph of geometry and flexibility

    def append_residue(self, residue):
        residue.protein = self
        residue.label = len(self.residues)
        self.residues.append(residue)
        
    def append_segment(self, segment):
        self.segments.append(segment)
        
    def make_graph(self):
        self.graph.parse(self.residues, self.all_flexible, self.segments)

        for d, indices in enumerate(self.graph.rigid_domains):
            for i in indices:
                self.residues[i].rigid_domain = d
                
        for b in self.graph.segment_bonds:
           self.residues[b.i].segment_bond = b.j
           self.residues[b.j].segment_bond = b.i
            


class Potential(object):
    def __init__(self, proteins, charge, radius, pair_potentials, torsions):
        self.proteins = proteins
        self.charge = charge
        self.radius = radius
        self.pair_potentials = pair_potentials
        self.torsions = torsions
        self.components = defaultdict(float)

    @staticmethod
    def distance(pos1, pos2, bounding_value=None):
        diffs = [c1 - c2 for (c1, c2) in zip(pos1, pos2)]
        if bounding_value is not None:
            diffs = [d - bounding_value*round(d/float(bounding_value)) for d in diffs]
        return sum(d**2 for d in diffs)**0.5

    @classmethod
    def from_filelist(cls, args):

        charge = {}
        radius = {}
        pair_potentials = defaultdict(dict)
        torsions = defaultdict(dict)
        proteins = []

        with open(args.aminoacids, 'r') as aa_file:
            acids = AMINOACID.findall(aa_file.read())

        for a, r, c in acids:
            amino_acid = a.upper()
            radius[amino_acid] = float(r)
            charge[amino_acid] = float(c)

        with open(args.potential, 'r') as pp_file:
            lines = [l for l in pp_file if not COMMENT.search(l)]

        amino_acids = [a.upper() for a in lines[0].strip().split()]

        for l in lines[1:]:
            chunks = l.strip().split()
            acid_i = chunks[0].upper()
            other_acids = amino_acids[amino_acids.index(acid_i):]
            potentials = [float(p) for p in chunks[1:]]

            for acid_j, potential_value in zip(other_acids, potentials):
                pair_potentials[acid_i][acid_j] = potential_value
                pair_potentials[acid_j][acid_i] = potential_value

        with open(args.torsions, 'r') as t_file:
            lines = [l for l in t_file if not COMMENT.search(l)]

        for l in lines:
            aa1, aa2, v, n, sigma = TORSION.match(l).groups()
            if aa2 not in torsions[aa1]:
                torsions[aa1][aa2] = {}
            torsions[aa1][aa2][int(n)] = {"V": float(v), "sigma": np.radians(float(sigma))}

        for filename in args.files:
            with open(filename) as pdbfile:
                protein = None

                for line in pdbfile:
                    ca_atom = ATOM.match(line)
                    if ca_atom:
                        amino_acid, chain, x, y, z = ca_atom.groups()
                        if not chain:
                            chain = 'A'

                        if not protein:
                            protein = Protein()
                        protein.append_residue(Residue(amino_acid, float(x), float(y), float(z), chain))

                    elif line.startswith("TER"):
                        if args.split_chains:
                            proteins.append(protein)
                            protein = None

                if protein:
                    proteins.append(protein)

        for mi in args.all_flexible:
            proteins[mi].all_flexible = True
            
        for s in args.segments:
            m, r = s.split(':')
            m = int(m)
            r = [int(i) for i in r.split(',')]
            proteins[m].append_segment(r)
            
        for p in proteins:
            p.make_graph()

        return cls(proteins, charge, radius, pair_potentials, torsions)

    def calculate(self, args):
        def debug(message):
            if args.debug_potential:
                print message
        residues = [r for p in self.proteins for r in p.residues]

        debug("UNBONDED")

        for i, ri in enumerate(residues):
            for rj in residues[i + 1:]:
                if ri.protein == rj.protein:
                    # ignore this pair if we're not calculating internal E at all or this specific molecule is rigid
                    if not args.internal or not ri.protein.segments:
                        continue

                    # ignore this pair if residues are close to each other on the backbone of the same chain
                    if ri.chain == rj.chain and abs(ri.label - rj.label) < 4:
                        continue
                    
                    # ignore this pair if residues are in the same rigid domain
                    if ri.rigid_domain is not None and ri.rigid_domain == rj.rigid_domain:
                        continue
                    
                    # ignore this pair if residues are bonded to each other
                    if ri.segment_bond == rj.label:
                        continue
                    
                    # ignore this pair if residues are close to each other because of a bond
                    if (ri.label, rj.label) in p.graph.indirect_neighbours:
                        continue

                r = self.distance(rj.position, ri.position, args.boundary) + EPS
                Eij = LAMBDA * (self.pair_potentials[ri.amino_acid][rj.amino_acid] - E0)
                sigmaij = (self.radius[ri.amino_acid] + self.radius[rj.amino_acid])/2.0
                r0ij = 2**(1/6.0) * sigmaij

                if Eij < 0:
                    LJ = 4 * np.abs(Eij) * ((sigmaij / r)**12 - (sigmaij / r)**6)
                elif r < r0ij:
                    LJ = 4 * Eij * ((sigmaij / r)**12 - (sigmaij / r)**6) + 2 * Eij
                else: # r >= r0ij
                    LJ = -4 * Eij * ((sigmaij / r)**12 - (sigmaij / r)**6)

                self.components["LJ"] += LJ * RT_TO_KCALMOL

                DH = self.charge[ri.amino_acid] * self.charge[rj.amino_acid] * np.exp(-r / XI) / r
                self.components["DH"] += DH * DH_CONVERSION_FACTOR
                
                debug("i: %d\tj: %d\tLJ: %f\tDH: %f" % (ri.label, rj.label, LJ * RT_TO_KCALMOL, DH * DH_CONVERSION_FACTOR))

        if args.internal:
            debug("BONDED")

            for p in self.proteins:

                for b in p.graph.bonds:
                    r1 = p.residues[b.i]
                    r2 = p.residues[b.j]
                    
                    r = np.linalg.norm(r1.position - r2.position)
                    
                    debug(r)
                    
                    bond_e =  K_SPRING * (r - R0)**2 / 2.0

                    self.components["bond"] += bond_e
                    
                    debug("bond %d-%d: %g" % (b.i, b.j, bond_e))
                    
                for a in p.graph.angles:
                    r1 = p.residues[a.i]
                    r2 = p.residues[a.j]
                    r3 = p.residues[a.k]
                    
                    
                    ba = r1.position - r2.position
                    bc = r3.position - r2.position
                    theta = np.arccos(np.dot(ba, bc) / np.linalg.norm(ba) / np.linalg.norm(bc))
                    print "THETA: %f" % theta
                    
                    angle_e = np.log(np.exp(-GAMMA_ANGLE * (K_ALPHA * (theta - THETA_ALPHA)**2 + EPSILON_ALPHA)) + np.exp(-GAMMA_ANGLE * K_BETA * (theta - THETA_BETA)**2)) / - GAMMA_ANGLE
                    print "ANGLE POTENTIAL: %f" % angle_e

                    self.components["angle"] += angle_e
                    
                    debug("angle %d-%d-%d: %g" % (a.i, a.j, a.k, angle_e))
                    
                for t in p.graph.torsions:
                    r1 = p.residues[t.i]
                    r2 = p.residues[t.j]
                    r3 = p.residues[t.k]
                    r4 = p.residues[t.l]
                    
                    b1 = r2.position - r1.position
                    b2 = r3.position - r2.position
                    b3 = r4.position - r3.position
                    b2xb3 = np.cross(b2, b3)
                    phi = np.arctan2((np.linalg.norm(b2) * np.dot(b1, b2xb3)), np.dot(np.cross(b1, b2), b2xb3))

                    aa1 = r2.amino_acid
                    aa2 = r3.amino_acid
                    
                    torsion_e = sum((1 + np.cos(n * phi - self.torsions[aa1][aa2][n]["sigma"])) * self.torsions[aa1][aa2][n]["V"] for n in range(1, 5))

                    self.components["torsion"] += torsion_e
                    
                    debug("torsion %d-%d-%d-%d: %g" % (t.i, t.j, t.k, t.l, torsion_e))

        self.components["total"] = 0
        self.components["total"] = sum(self.components.values())

    def __str__(self):
        return ", ".join("%s: %f" % (k, self.components[k]) for k in ("total", "LJ", "DH", "bond", "angle", "torsion"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("aminoacids", help="Amino acid data file")
    parser.add_argument("potential", help="Pair potential file")
    parser.add_argument("torsions", help="Torsion data file")
    parser.add_argument("files", help="Input PDB files", nargs="+")
    parser.add_argument("-s" "--segment", help="Flexible segment (format: 'mi:ri,rj,rk,(...),rz')", action="append", dest="segments", default=[])
    parser.add_argument("-f" "--all-flexible", help="Molecule is all flexible", type=int, action="append", dest="all_flexible", default=[])
    parser.add_argument("-i", "--internal", help="Calculate internal molecule potential", action="store_true")
    parser.add_argument("-c", "--split-chains", help="Treat each chain as a separate protein, even if they are in the same PDB file", action="store_true")
    parser.add_argument("-d", "--debug-potential", help="Print details about the potential calculation", action="store_true")
    parser.add_argument("-b", "--boundary", help="Periodic boundary", type=float, default=None)
    parser.add_argument("-l", "--lj", help="Type of LJ potential", type=str, default="normal", choices=['normal', 'repulsive', 'off'])
    args = parser.parse_args()
    
    if args.lj == "repulsive":
        E0 = 0.0001
    elif args.lj == "off":
        LAMBDA = 0

    potential = Potential.from_filelist(args)
    potential.calculate(args)
    print potential
