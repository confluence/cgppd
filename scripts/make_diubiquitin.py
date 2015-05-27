#!/usr/bin/env python

import sys
import numpy as np
import math
import re
import argparse
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

    ATOM_FORMAT = "ATOM  {:>5}  CA  {} {}{:>4}    {:8.3f}{:8.3f}{:8.3f}  0.00  0.00\n" # atom id, amino acid, chain id, residue id, x, y, z
    TER_FORMAT = "TER   {:>5}      {} {}{:>4}\n"

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
            lines.append(self.ATOM_FORMAT.format(offset + i + 1, r.amino_acid, chain_id, i + 1, r.pos[0], r.pos[1], r.pos[2]))

        i = len(self.residues)
        r = self.residues[-1]
        lines.append(self.TER_FORMAT.format(offset + i, r.amino_acid, chain_id, i))

        return lines

    @property
    def centre(self):
        return sum(r.pos for r in self.residues)/float(len(self.residues))

    @property
    def end(self):
        return self.residues[-1].pos

    
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

        lines.append("END\n")

        return lines


class Quaternion(object):
    def __init__(self, q):
        self.q = np.array([q[0], q[1], q[2], q[3]], dtype=np.float64)

    @classmethod
    def from_angle_axis(cls, angle, axis):
        q = cls.from_vector(axis).q
        l = np.linalg.norm(q)
        if l > np.finfo(float).eps * 4.0:
            q *= math.sin(angle/2.0) / l
        q[0] = math.cos(angle/2.0)
        return cls(q)

    @classmethod
    def from_vector(cls, v):
        vn = v / np.linalg.norm(v)
        return cls(np.array([0.0, vn[0], vn[1], vn[2]]))

    def conjugate(self):
        c = np.array(self.q, copy=True)
        np.negative(c[1:], c[1:])
        return Quaternion(c)

    def multiply(self, other):
        w0, x0, y0, z0 = self.q
        w1, x1, y1, z1 = other.q

        return Quaternion([
            -x1*x0 - y1*y0 - z1*z0 + w1*w0,
            x1*w0 + y1*z0 - z1*y0 + w1*x0,
            -x1*z0 + y1*w0 + z1*x0 + w1*y0,
            x1*y0 - y1*x0 + z1*w0 + w1*z0
        ])

    def rotate(self, v):
        return self.multiply(Quaternion.from_vector(v)).multiply(self.conjugate())[1:]


def test_quaternion_from_angle_axis():
    q = Quaternion.from_angle_axis(0.123, [1, 0, 0])
    assert np.allclose(q.q, [0.99810947, 0.06146124, 0, 0])

def test_quaternion_multiply():
    q = Quaternion([4, 1, -2, 3]).multiply(Quaternion([8, -5, 6, 7]))
    assert np.allclose(q.q, [28, 20, 30, 56]) # is this right?

    # TODO a rotation test, and then replace the rotation from the visual library


UBIQUITIN_PDB = """ATOM     25  CA  MET     1      19.902   3.678  -0.899  0.00  0.00      B
ATOM     26  CA  GLN     2      17.074   3.470   1.639  0.00  0.00      B
ATOM     27  CA  ILE     3      13.397   4.446   1.531  0.00  0.00      B
ATOM     28  CA  PHE     4      10.485   4.314   3.980  0.00  0.00      B
ATOM     29  CA  VAL     5       7.237   2.625   2.988  0.00  0.00      B
ATOM     30  CA  LYS     6       4.171   4.223   4.562  0.00  0.00      B
ATOM     31  CA  THR     7       1.034   2.108   4.839  0.00  0.00      B
ATOM     32  CA  LEU     8      -2.335   3.794   4.306  0.00  0.00      B
ATOM     33  CA  THR     9      -3.231   2.679   7.835  0.00  0.00      B
ATOM     34  CA  GLY    10      -0.265   4.548   9.294  0.00  0.00      B
ATOM     35  CA  LYS    11       3.028   2.658   9.476  0.00  0.00      B
ATOM     36  CA  THR    12       6.434   3.592   8.055  0.00  0.00      B
ATOM     37  CA  ILE    13       9.111   0.946   7.468  0.00  0.00      B
ATOM     38  CA  THR    14      12.743   1.080   6.348  0.00  0.00      B
ATOM     39  CA  LEU    15      13.403  -0.504   2.956  0.00  0.00      B
ATOM     40  CA  GLU    16      16.797  -0.602   1.246  0.00  0.00      B
ATOM     41  CA  VAL    17      16.715  -0.687  -2.562  0.00  0.00      B
ATOM     42  CA  GLU    18      19.320  -0.022  -5.255  0.00  0.00      B
ATOM     43  CA  SER    19      18.584   2.969  -7.488  0.00  0.00      B
ATOM     44  CA  SER    20      19.324   0.871 -10.574  0.00  0.00      B
ATOM     45  CA  ASP    21      16.904  -1.835  -9.420  0.00  0.00      B
ATOM     46  CA  THR    22      13.387  -1.372 -10.796  0.00  0.00      B
ATOM     47  CA  ILE    23      10.175  -1.159  -8.785  0.00  0.00      B
ATOM     48  CA  ASP    24       9.505  -4.877  -9.267  0.00  0.00      B
ATOM     49  CA  ASN    25      12.525  -5.598  -7.067  0.00  0.00      B
ATOM     50  CA  VAL    26      11.218  -3.260  -4.368  0.00  0.00      B
ATOM     51  CA  LYS    27       7.824  -4.974  -4.378  0.00  0.00      B
ATOM     52  CA  SER    28       9.572  -8.333  -4.052  0.00  0.00      B
ATOM     53  CA  LYS    29      11.341  -7.233  -0.869  0.00  0.00      B
ATOM     54  CA  ILE    30       8.106  -5.750   0.471  0.00  0.00      B
ATOM     55  CA  GLN    31       6.369  -9.072  -0.146  0.00  0.00      B
ATOM     56  CA  ASP    32       9.001 -11.077   1.714  0.00  0.00      B
ATOM     57  CA  LYS    33       9.002  -8.736   4.725  0.00  0.00      B
ATOM     58  CA  GLU    34       5.381  -7.567   5.016  0.00  0.00      B
ATOM     59  CA  GLY    35       3.801 -10.435   3.089  0.00  0.00      B
ATOM     60  CA  ILE    36       1.752  -8.691   0.409  0.00  0.00      B
ATOM     61  CA  PRO    37       1.482  -9.402  -3.368  0.00  0.00      B
ATOM     62  CA  PRO    38       3.201  -6.805  -5.633  0.00  0.00      B
ATOM     63  CA  ASP    39       0.159  -6.442  -7.894  0.00  0.00      B
ATOM     64  CA  GLN    40      -2.135  -5.909  -4.904  0.00  0.00      B
ATOM     65  CA  GLN    41      -0.003  -3.205  -3.278  0.00  0.00      B
ATOM     66  CA  ARG    42       0.488   0.075  -5.144  0.00  0.00      B
ATOM     67  CA  LEU    43       3.020   2.771  -4.253  0.00  0.00      B
ATOM     68  CA  ILE    44       2.387   6.467  -4.892  0.00  0.00      B
ATOM     69  CA  PHE    45       4.813   9.397  -4.968  0.00  0.00      B
ATOM     70  CA  ALA    46       5.090  12.849  -6.545  0.00  0.00      B
ATOM     71  CA  GLY    47       1.467  12.546  -7.608  0.00  0.00      B
ATOM     72  CA  LYS    48       2.314   9.475  -9.686  0.00  0.00      B
ATOM     73  CA  GLN    49       1.875   5.743  -9.099  0.00  0.00      B
ATOM     74  CA  LEU    50       5.141   3.802  -9.350  0.00  0.00      B
ATOM     75  CA  GLU    51       5.076   1.333 -12.247  0.00  0.00      B
ATOM     76  CA  ASP    52       6.348  -2.212 -11.669  0.00  0.00      B
ATOM     77  CA  GLY    53       8.857  -1.828 -14.497  0.00  0.00      B
ATOM     78  CA  ARG    54      10.850   1.389 -14.173  0.00  0.00      B
ATOM     79  CA  THR    55      14.217   2.226 -12.605  0.00  0.00      B
ATOM     80  CA  LEU    56      14.234   4.507  -9.551  0.00  0.00      B
ATOM     81  CA  SER    57      16.912   6.658 -11.203  0.00  0.00      B
ATOM     82  CA  ASP    58      14.397   8.052 -13.699  0.00  0.00      B
ATOM     83  CA  TYR    59      12.055   8.936 -10.838  0.00  0.00      B
ATOM     84  CA  ASN    60      14.947  10.696  -9.093  0.00  0.00      B
ATOM     85  CA  ILE    61      13.548   9.712  -5.701  0.00  0.00      B
ATOM     86  CA  GLN    62      15.069  11.590  -2.765  0.00  0.00      B
ATOM     87  CA  LYS    63      17.110   9.644  -0.223  0.00  0.00      B
ATOM     88  CA  GLU    64      14.968   8.133   2.544  0.00  0.00      B
ATOM     89  CA  SER    65      11.745   9.187   0.804  0.00  0.00      B
ATOM     90  CA  THR    66       8.317   8.042   1.987  0.00  0.00      B
ATOM     91  CA  LEU    67       6.122   6.040  -0.403  0.00  0.00      B
ATOM     92  CA  HIS    68       2.329   5.919  -0.046  0.00  0.00      B
ATOM     93  CA  LEU    69       0.877   2.402   0.094  0.00  0.00      B
ATOM     94  CA  VAL    70      -2.444   1.677  -1.645  0.00  0.00      B
ATOM     95  CA  LEU    71      -4.024  -1.776  -1.826  0.00  0.00      B
ATOM     96  CA  ARG    72      -6.216  -2.372  -4.877  0.00  0.00      B
ATOM     97  CA  LEU    73      -9.833  -3.079  -3.935  1.00  0.00      B
ATOM     98  CA  ARG    74     -11.071  -6.607  -4.655  1.00  0.00      B
ATOM     99  CA  GLY    75     -13.233  -6.905  -7.755  1.00  0.00      B
ATOM    100  CA  GLY    76     -15.901  -8.793  -5.845  1.00  0.00      B"""


def make_diubiquitin(binding_resid, fake_sidechain=False):
    ubiquitin = Chain.from_pdb(UBIQUITIN_PDB.split('\n'))
    res = ubiquitin.residues[binding_resid - 1]
    sitename = "%s_%d" % (res.amino_acid.lower(), binding_resid)
    site = res.pos

    # target vector = binding site -> centre

    v_t = ubiquitin.centre - site

    # starting vector = centre -> end of tail

    v_s = ubiquitin.end - ubiquitin.centre

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

    tail_offset = 3.8 * (trans/np.linalg.norm(trans))

    if fake_sidechain:
        tail_offset *= 2 # make room for the fake alanine

    trans += tail_offset

    # create translated and rotated ubiquitin

    second_ubiquitin = Chain.copy(ubiquitin)
    second_ubiquitin.rotate(theta, axis, ubiquitin.centre)
    second_ubiquitin.translate(trans)

    # add the alanine

    if fake_sidechain:
        alanine = Residue("ALA", site + 0.5 * tail_offset)
        second_ubiquitin.residues.append(alanine)

    # output PDB

    diubiquitin = Molecule([second_ubiquitin, ubiquitin])
    if fake_sidechain:
        sitename += "_sidechain"
    output_filename = "diubiquitin_%s.pdb" % sitename

    return diubiquitin, output_filename


def test_make_diubiquitin():
    resids = (48, 63, 1)
    results = [make_diubiquitin(resid, False) for resid in resids]
    molecules, filenames = zip(*results)
    
    assert filenames == ("diubiquitin_lys_48.pdb", "diubiquitin_lys_63.pdb", "diubiquitin_met_1.pdb")

    for resid, molecule in zip(resids, molecules):
        assert len(molecule.chains) == 2

        ubq1, ubq2 = molecule.chains
        assert len(ubq1.residues) == 76
        assert len(ubq2.residues) == 76

        bond_length = np.linalg.norm(ubq1.residues[76 - 1].pos - ubq2.residues[resid - 1].pos)
        assert abs(bond_length - 3.8) < 0.00001
        

def test_make_diubiquitin_with_fake_sidechain():
    resids = (48, 63, 1)
    results = [make_diubiquitin(resid, True) for resid in resids]
    molecules, filenames = zip(*results)
    
    assert filenames == ("diubiquitin_lys_48_sidechain.pdb", "diubiquitin_lys_63_sidechain.pdb", "diubiquitin_met_1_sidechain.pdb")

    for resid, molecule in zip(resids, molecules):
        assert len(molecule.chains) == 2

        ubq1, ubq2 = molecule.chains
        assert len(ubq1.residues) == 77
        assert len(ubq2.residues) == 76

        bond_length = np.linalg.norm(ubq1.residues[76 - 1].pos - ubq2.residues[resid - 1].pos)
        tail_ala_length = np.linalg.norm(ubq1.residues[76 - 1].pos - ubq1.residues[77 - 1].pos)
        ala_site_length = np.linalg.norm(ubq1.residues[77 - 1].pos - ubq2.residues[resid - 1].pos)
        
        assert abs(bond_length - 3.8 * 2) < 0.00001
        assert abs(tail_ala_length - 3.8) < 0.00001
        assert abs(ala_site_length - 3.8) < 0.00001

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate diubiquitin PDBs")
    parser.add_argument("-f", "--fake-sidechain", help="Create a fake side-chain by inserting an alanine between the bonded residues", action="store_true")
    args = parser.parse_args()

    for resid in (48, 63, 1):
        diubiquitin, output_filename = make_diubiquitin(resid, args.fake_sidechain)
        with open(output_filename, 'w') as outputfile:
            outputfile.writelines(diubiquitin.to_pdb())

# TODO: check if we're actually doing quaternion rotation!
