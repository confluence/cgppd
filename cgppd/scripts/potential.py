#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
from collections import defaultdict
import numpy as np
import argparse

COMMENT = re.compile(r"^\s*#")
ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) [A-Z]? *\d+ *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(\d+\.\d{2}).*')


class Residue(object):
    def __init__(self, amino_acid, x, y, z, flexible, molecule_no):
        self.amino_acid = amino_acid
        self.flexible = flexible
        self.position = np.array((x, y, z))
        self.molecule_no = molecule_no


class Protein(object):
    def __init__(self):
        self.residues = []

    def append_residue(self, residue):
        self.residues.append(residue)


class Potential(object):

    def __init__(self, proteins, pair_potentials):
        self.proteins = proteins
        self.pair_potentials = pair_potentials
        self.components = {}

    @classmethod
    def from_filelist(cls, potentialfilename, *pdbfilenames):
        pair_potentials = defaultdict(dict)
        proteins = []

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
            print "file %s" % filename
            with open(filename) as pdbfile:
                protein = None

                for line in pdbfile:
                    ca_atom = ATOM.match(line)
                    if ca_atom:
                        amino_acid, x, y, z, flexible = ca_atom.groups()

                        if not protein:
                            protein = Protein()
                        protein.append_residue(Residue(amino_acid, float(x), float(y), float(z), bool(float(flexible)), len(proteins)))

                    elif line.startswith("TER"):
                        proteins.append(protein)
                        protein = None

                if protein:
                    proteins.append(protein)

        return cls(proteins, pair_potentials)

    @property
    def residues(self):
        return [r for p in self.proteins for r in p.residues]


# LAMBDA
# E0

    @property
    def LJ(self, include_internal):

        if "LJ" not in self.components:

            for i, ri in enumerate(self.residues):
                for rj in self.residues[i + 1:]:
                    if not include_internal and ri.molecule_no == rj.molecule_no:
                        continue

                    r = np.linalg.norm(rj.position - ri.position)

    #double Eij(lambda * (AminoAcidsData.LJpotentials[ri.aminoAcidIndex][rj.aminoAcidIndex] - e0));

    #// sigmaij = (sigmai + sigmaj) / 2
    #double sigmaij(0.5f * (ri.vanderWaalRadius + rj.vanderWaalRadius));
    #double LJtmp(powf(sigmaij / r, 6.0f)); // (sigmaij/r)^6

    #// if attractive interactions (Eij < 0), or repulsive interactions (Eij > 0) and r >= r0ij:
    #// uij(r) = -4Eij( (sigmaij/r)^12 - (sigmaij/r)^6 )
    #double LJ(-4.0f * Eij * LJtmp * (LJtmp - 1.0f));

    #// if repulsive interactions (Eij > 0) and r < r0ij:
    #// uij(r) = 4Eij( (sigmaij/r)^12 - (sigmaij/r)^6 ) + 2Eij
    #if (Eij > 0.0f && r < sigmaij * r0_constant)
    #{
        #LJ = -LJ + 2.0f * Eij;
    #}

        return self.components["LJ"]

    @property
    def DH(self):
        pass # TODO

    @property
    def bond(self):
        pass # TODO

    @property
    def angle(self):
        pass # TODO

    @property
    def torsion(self):
        pass # TODO

    @property
    def total(self):
        pass # TODO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation output from cgppd")
    parser.add_argument("potential", help="Pair potential file")
    parser.add_argument("files", help="Input PDB files", nargs="+")
    args = parser.parse_args()

    potential = Potential.from_filelist(args.potential, *args.files)
