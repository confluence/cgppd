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

binding_sites = {
    "LYS 48": ubiquitin[48 - 1],
    "LYS 63": ubiquitin[63 - 1],
    "MET 1": ubiquitin[1 - 1],
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

    print "%s: translation %s rotation axis %s theta %s" % (sitename, trans, axis, theta)
