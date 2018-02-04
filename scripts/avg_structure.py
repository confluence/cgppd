#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import os
import glob
import numpy as np

clusterfilename = sys.argv[1]
dirname = os.path.dirname(clusterfilename)

def sampleglob():
    return glob.iglob(os.path.join(dirname, "pdb", "sample*_303.8K_*.pdb"))

ATOM = re.compile(r'ATOM *(\d+) *CA *[A-Z]{3} *[A-Z] *\d+ *(-?\d+.\d{3}) *(-?\d+.\d{3}) *(-?\d+.\d{3}) *0.00 *0.00')

# get cluster frames

clusters = []

with open(clusterfilename) as clusterfile:
    for line in clusterfile:
        clusters.append([int(f) for f in line.split()])
        
print("Found %d clusters." % len(clusters))

reverse_mapping = {f:i for (i, c) in enumerate(clusters) for f in c}

# calculate centres

centres = [np.zeros((153, 3)) for c in clusters]

for sample_id, samplefilename in enumerate(sampleglob()):
    cluster_id = reverse_mapping.get(sample_id + 1)
    if cluster_id is None:
        continue
    
    num_samples = len(clusters[cluster_id])
    
    with open(samplefilename) as samplefile:
        for line in samplefile:
            m = ATOM.search(line)
            if m:
                r, x, y, z = m.groups()
                r, x, y, z = int(r), float(x), float(y), float(z)
                centres[cluster_id][r - 1] += np.array((x, y, z)) / num_samples
        
print("Found centres.")

# calculate minumum distance

rmsds = [[] for c in clusters]

for sample_id, samplefilename in enumerate(sampleglob()):
    cluster_id = reverse_mapping.get(sample_id + 1)
    if cluster_id is None:
        continue
    
    sq_atom_distances = []
    
    with open(samplefilename) as samplefile:
        for line in samplefile:
            m = ATOM.search(line)
            if m:
                r, x, y, z = m.groups()
                r, x, y, z = int(r), float(x), float(y), float(z)
                sq_atom_distances.append(np.linalg.norm(np.array((x, y, z)) - centres[cluster_id][r - 1])**2)
                
    rmsd = np.sqrt(sum(sq_atom_distances)/153)
    rmsds[cluster_id].append((rmsd, sample_id))
    

print("Frames closest to centre, per cluster:")
for distances in rmsds:
    print([f for d, f in sorted(distances)[:10]])
