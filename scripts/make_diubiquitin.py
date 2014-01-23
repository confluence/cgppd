#!/usr/bin/env python

import sys
import numpy as np
import re

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
print "centre %s" % centre
end = ubiquitin[76 - 1]
print "end %s" % end
lys_48 = ubiquitin[48 - 1]
print "lys 48 %s" % lys_48

# target vector = binding site -> centre

v_t = centre - lys_48

print "v_t %s" % v_t

# starting vector = centre -> end of tail

v_s = end - centre

print "v_s %s" % v_s

# calculate rotation axis and angle from these vectors

A = np.cross(v_s, v_t)
axis = A/np.linalg.norm(A)

sin_theta = np.linalg.norm(A) / (np.linalg.norm(v_s) * np.linalg.norm(v_t))
theta = np.arcsin(sin_theta)

print "axis %s theta %s" % (axis, theta)

# apply rotation to end of tail (relative to centre; i.e. v_s)

a = np.cos(theta/2)
b,c,d = -axis * np.sin(theta/2)
rot_matrix = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])



rot_v_s = np.dot(rot_matrix, v_s)
print "rot_v_s %s" % rot_v_s

# calculate translation vector: end of rotated tail -> binding site (-v_t!)

trans = -v_t - rot_v_s

print "trans %s" % trans

# add padding of 3.8 A in the same direction

padding = 3.8 * (trans/np.linalg.norm(trans))

print "padding %s" % padding

trans += padding

print "trans %s" % trans

# output translation and rotation

print "LYS 48: translation %s rotation axis %s theta %s" % (trans, axis, theta)
