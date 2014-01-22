#!/usr/bin/env python

import re
import sys
import numpy as np

ATOM = re.compile('ATOM *\d+ *CA *([A-Z]{3}) ([A-Z ]) *(\d+) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *(-?\d+\.\d{3}) *\d+\.\d{2} *\d+\.\d{2}')

with open(sys.argv[1]) as pdbfile:
    last_pos = None

    for line in pdbfile:
        m = ATOM.match(line)
        if not m:
            continue

        amino_acid, chain_id, index, x, y, z = m.groups()
        pos = np.array((float(x), float(y), float(z)))

        if last_pos is not None:
            print np.linalg.norm(pos - last_pos)

        last_pos = pos
