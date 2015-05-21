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



"""
ATOM     25  CA  MET     1      19.902   3.678  -0.899  0.00  0.00      B
ATOM     26  CA  GLN     2      17.074   3.470   1.639  0.00  0.00      B
ATOM     27  CA  ILE     3      13.397   4.446   1.531  0.00  0.00      B
ATOM     28  CA  PHE     4      10.485   4.314   3.980  0.00  0.00      B
ATOM     29  CA  VAL     5       7.237   2.625   2.988  0.00  0.00      B
ATOM     30  CA  LYS     6       4.171   4.223   4.562  0.00  0.00      B
"""

"""
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34      A1   C
ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65      A1   O
ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88      A1   C
ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41      A1   C
ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64      A1   C
ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11      A1   C
ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58      A1   C
ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25      A1   C
"""
        

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

    def rotate(self, theta, axis):
        for r in self.residues:
            r.pos = np.array(vector(r.pos).rotate(theta, axis))

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


        

#ubiquitin = []

#with open(sys.argv[1]) as pdbfile:
    #for line in pdbfile:
        #m = ATOM.match(line)

        #if not m:
            #continue

        #amino_acid, chain_id, index, x, y, z = m.groups()
        #ubiquitin.append(np.array((float(x), float(y), float(z))))

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

    # output translation and rotation

    #print "# %s:\nt(%f, %f, %f) r(%f, %f, %f, %f) %s" % (sitename, trans[0], trans[1], trans[2], axis[0], axis[1], axis[2], theta, sys.argv[1])

    # create translated and rotated ubiquitin

    second_ubiquitin = Chain.copy(ubiquitin)
    second_ubiquitin.rotate(theta, axis)
    second_ubiquitin.translate(trans)

    # TODO this is where we will add the alanine

    # output PDB

    molecule = Molecule([second_ubiquitin, ubiquitin])
    output_filename = "diubiquitin_%s.pdb" % sitename

    with open(output_filename, 'w') as outputfile:
        outputfile.writelines(molecule.to_pdb())


# TODO: optionally add fake alanine sidechain to end of tail
#       pad the bond with 2 * 3.8 instead and put the alanine in the middle
# TODO: actually apply the rotation to the ubiquitin and output a PDB directly
