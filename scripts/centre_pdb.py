#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import numpy as np

if len(sys.argv) < 2:
    sys.exit("Please enter a file parameter.")

filename = sys.argv[1]

# First pass to calculate geometric centre

total = np.array([0.0, 0.0, 0.0])
num_atoms = 0

with open(filename, "r") as pdbfile:
    for line in pdbfile:
        if line.startswith("ATOM"):
            pos = np.array([float(n) for n in re.findall('(-?\d+\.\d{3})', line[30:-12])])
            total += pos
            num_atoms += 1
            

centre = total / num_atoms

with open(filename, "r") as pdbfile:
    for line in pdbfile:
        if line.startswith("ATOM"):
            pos = np.array([float(n) for n in re.findall('(-?\d+\.\d{3})', line[30:-12])])
            pos -= centre
            new_line = line[:30] + ("%8.3f%8.3f%8.3f" % (pos[0], pos[1], pos[2])) + line[-12:]
            sys.stdout.write(new_line)
        else:
            sys.stdout.write(line)
