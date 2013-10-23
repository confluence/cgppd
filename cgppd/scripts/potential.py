#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
from collections import defaultdict
import numpy as np
import argparse

E0 = -2.27
LAMBDA = 0.159
XI = 10.0
LJ_CONVERSION_FACTOR = 0.3507221006079
BTU_to_J = 0.948
E_CHARGE = 1.602176487 # what about the e-19?
KBT_TO_KCALMOL = 1.0 / (294.0 * 8.314472 / 4184.0)

COMMENT = re.compile(r"^\s*#")
ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) [A-Z]? *\d+ *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(\d+\.\d{2}).*')
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

class Residue(object):
    def __init__(self, amino_acid, x, y, z, flexible, molecule_no, residue_no):
        self.amino_acid = amino_acid
        self.flexible = flexible
        self.position = np.array((x, y, z))
        self.molecule_no = molecule_no
        self.residue_no = residue_no


class Protein(object):
    def __init__(self):
        self.residues = []
        self.segments = [] # flexible segments

    def append_residue(self, residue):
        self.residues.append(residue)

    def find_segments(self):
        segment = None
        for r in self.residues:
            if r.flexible:
                if not segment:
                    segment = []
                segment.append(r)

            if not r.flexible:
                if segment:
                    self.segments.append(segment)
                    segment = None


class Potential(object):
    def __init__(self, proteins, charge, radius, pair_potentials):
        self.proteins = proteins
        self.charge = charge
        self.radius = radius
        self.pair_potentials = pair_potentials
        self.components = defaultdict(float)

    @classmethod
    def from_filelist(cls, aminoacidfilename, potentialfilename, *pdbfilenames):
        charge = {}
        radius = {}
        pair_potentials = defaultdict(dict)
        proteins = []

        with open(aminoacidfilename, 'r') as aa_file:
            acids = AMINOACID.findall(aa_file.read())

        for a, r, c in acids:
            amino_acid = a.upper()
            radius[amino_acid] = float(r)
            charge[amino_acid] = float(c)

        with open(potentialfilename, 'r') as pp_file:
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

        for filename in pdbfilenames:
            with open(filename) as pdbfile:
                protein = None
                residue_no = 0

                for line in pdbfile:
                    ca_atom = ATOM.match(line)
                    if ca_atom:
                        amino_acid, x, y, z, flexible = ca_atom.groups()

                        if not protein:
                            protein = Protein()
                        protein.append_residue(Residue(amino_acid, float(x), float(y), float(z), bool(float(flexible)), len(proteins), residue_no))
                        residue_no += 1

                    elif line.startswith("TER"):
                        proteins.append(protein)
                        protein = None
                        residue_no = 0

                if protein:
                    proteins.append(protein)

        for protein in proteins:
            protein.find_segments()

        return cls(proteins, charge, radius, pair_potentials)

    def calculate(self, include_internal):
        residues = [r for p in self.proteins for r in p.residues]

        for i, ri in enumerate(residues):
            for rj in residues[i + 1:]:
                if not include_internal and ri.molecule_no == rj.molecule_no:
                    continue

                if ri.molecule_no == rj.molecule_no and np.abs(ri.residue_no - rj.residue_no) < 4:
                    continue

                r = np.linalg.norm(rj.position - ri.position)
                Eij = LAMBDA * (self.pair_potentials[ri.amino_acid][rj.amino_acid] - E0)
                sigmaij = (self.radius[ri.amino_acid] + self.radius[rj.amino_acid])/2.0
                r0ij = 2**(1/6) * sigmaij

                if Eij < 0:
                    LJ = 4 * np.abs(Eij) * ((sigmaij / r)**12 - (sigmaij / r)**6)
                elif r < r0ij:
                    LJ = 4 * Eij * ((sigmaij / r)**12 - (sigmaij / r)**6) + 2 * Eij
                else: # r >= r0ij
                    LJ = -4 * Eij * ((sigmaij / r)**12 - (sigmaij / r)**6)

                self.components["LJ"] += LJ * LJ_CONVERSION_FACTOR * KBT_TO_KCALMOL

                DH = self.charge[ri.amino_acid] * self.charge[rj.amino_acid] * np.exp(-r / XI) / r
                self.components["DH"] += DH * BTU_to_J * E_CHARGE * E_CHARGE * KBT_TO_KCALMOL

                # TODO: why are the unit conversions what they are?
                # TODO: WTF happened to the 4piD?



        self.components["total"] = 0
        self.components["total"] = sum(self.components.values())

    def __str__(self):
        return "\n".join("%s: %g" % (k, self.components[k]) for k in ("LJ", "DH", "bond", "angle", "torsion", "total"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("aminoacids", help="Amino acid data file")
    parser.add_argument("potential", help="Pair potential file")
    parser.add_argument("files", help="Input PDB files", nargs="+")
    parser.add_argument("-i", "--internal", help="Calculate internal molecule potential", action="store_true")
    args = parser.parse_args()

    potential = Potential.from_filelist(args.aminoacids, args.potential, *args.files)
    potential.calculate(args.internal)
    print potential
