#!/usr/bin/env python

import sys
import numpy as np
import re
from visual import vector


class Residue(object):
    def __init__(self, amino_acid, pos):
        self.amino_acid = amino_acid
        self.pos = pos

    @classmethod
    def copy(cls, residue):
        return cls(residue.amino_acid, np.copy(residue.pos))


class Chain(object):
    ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z ]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')

    FORMAT = "ATOM  {:>5}  CA  {} {}{:>4}    {:8.3f}{:8.3f}{:8.3f}  0.00  0.00\n" # atom id, amino acid, chain id, residue id, x, y, z

    def __init__(self, residues):
        self.residues = residues

    @classmethod
    def from_pdb(cls, pdbfile):
        residues = []
        
        for line in pdbfile:
            m = cls.ATOM.match(line)

            if not m:
                continue

            amino_acid, chain_id, index, x, y, z = m.groups()
            residues.append(Residue(amino_acid, np.array((float(x), float(y), float(z)))))

        return cls(residues)

    @classmethod
    def copy(cls, chain):
        return cls([Residue.copy(r) for r in chain.residues])

    def rotate(self, theta, axis, centre):
        for r in self.residues:
            rel_pos = r.pos - centre
            rel_pos = np.array(vector(rel_pos).rotate(theta, axis))
            r.pos = rel_pos + centre

    def translate(self, translation_vector):
        for r in self.residues:
            r.pos += translation_vector

    def to_pdb(self, chain_id, offset):
        lines = []

        for i, r in enumerate(self.residues):
            lines.append(self.FORMAT.format(offset + i + 1, r.amino_acid, chain_id, i + 1, r.pos[0], r.pos[1], r.pos[2]))

        return lines

    
class Molecule(object):
    def __init__(self, chains):
        self.chains = chains

    def to_pdb(self):
        offset = 0
        chain_id = 'A'
        lines = []
        
        for c in self.chains:
            lines.extend(c.to_pdb(chain_id, offset))
            offset += len(c.residues)
            chain_id = chr(ord(chain_id) + 1)

        lines.append("TER\n")

        return lines


with open(sys.argv[1]) as pdbfile:
    ubiquitin = Chain.from_pdb(pdbfile)


centre = sum(r.pos for r in ubiquitin.residues)/float(len(ubiquitin.residues))
end = ubiquitin.residues[76 - 1].pos

binding_sites = {
    "lys_48": ubiquitin.residues[48 - 1].pos,
    "lys_63": ubiquitin.residues[63 - 1].pos,
    "met_1": ubiquitin.residues[1 - 1].pos,
}

for sitename, site in binding_sites.iteritems():
    # target vector = binding site -> centre

    v_t = centre - site

    # starting vector = centre -> end of tail

    v_s = end - centre

    # calculate rotation axis and angle from these vectors

    A = np.cross(v_s, v_t)
    axis = A/np.linalg.norm(A)

    sin_theta = np.linalg.norm(A) / (np.linalg.norm(v_s) * np.linalg.norm(v_t))
    theta = np.arcsin(sin_theta)

    # apply rotation to end of tail (relative to centre; i.e. v_s)

    rot_v_s = np.array(vector(v_s).rotate(theta, axis))

    # calculate translation vector: end of rotated tail -> binding site (-v_t!)

    trans = -v_t - rot_v_s

    # add padding of 3.8 A in the same direction

    trans += 3.8 * (trans/np.linalg.norm(trans))

    # create translated and rotated ubiquitin

    second_ubiquitin = Chain.copy(ubiquitin)
    second_ubiquitin.rotate(theta, axis, centre)
    second_ubiquitin.translate(trans)

    # TODO this is where we will add the alanine

    # output PDB

    molecule = Molecule([second_ubiquitin, ubiquitin])
    output_filename = "diubiquitin_%s.pdb" % sitename

    with open(output_filename, 'w') as outputfile:
        outputfile.writelines(molecule.to_pdb())


# TODO: optionally add fake alanine sidechain to end of tail
#       pad the bond with 2 * 3.8 instead and put the alanine in the middle


# TODO: build in a test to measure the distance between the tail and the binding site.
