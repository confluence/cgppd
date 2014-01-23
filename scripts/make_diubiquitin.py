#!/usr/bin/env python

import sys
import numpy as np
import re
from visual import vector

ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z ]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')

ubiquitin = []

with open(sys.argv[1]) as pdbfile:
    for line in pdbfile:
        m = ATOM.match(line)

        if not m:
            continue

        amino_acid, chain_id, index, x, y, z = m.groups()
        ubiquitin.append(np.array((float(x), float(y), float(z))))

centre = sum(ubiquitin)/float(len(ubiquitin))
end = ubiquitin[76 - 1]
lys_48 = ubiquitin[48 - 1]

# target vector = binding site -> centre

v_t = centre - lys_48

# starting vector = centre -> end of tail

v_s = end - centre

# calculate rotation axis and angle from these vectors

A = np.cross(v_s, v_t)
axis = A/np.linalg.norm(A)

cos_theta = np.dot(v_s, v_t) / (np.linalg.norm(v_s) * np.linalg.norm(v_t))
theta = np.arccos(cos_theta)

# apply rotation to end of tail (relative to centre; i.e. v_s)

rot_v_s = np.array(vector(v_s).rotate(theta, axis))

# calculate translation vector: end of rotated tail -> binding site

trans = v_t - rot_v_s

# add padding of 8.8 A in the same direction

padding = 8.8 * trans/np.linalg.norm(trans)
trans += padding

# output translation and rotation

print "LYS 48: translation %s rotation axis %s angle %s" % (trans, axis, theta)
