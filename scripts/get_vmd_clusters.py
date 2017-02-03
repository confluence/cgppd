#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import os

NEW_MOLECULE = re.compile("mol new (.*\.pdb)")
CLUSTER = re.compile("mol drawframes top \d {(\d.*)}")

clusters = {}
current_cluster = None

with open(sys.argv[1], "r") as vis_file:
    for line in vis_file:
        nm = NEW_MOLECULE.match(line)
        if nm:
            current_cluster = nm.group(1).split('/')[-2]
            clusters[current_cluster] = []
            continue
            
        cl = CLUSTER.match(line)
        if cl:
            clusters[current_cluster].append(cl.group(1))

for dirname, frames in clusters.iteritems():
    if os.path.isdir(dirname):
        filename = "%s/clusters.txt" % dirname
    else:
        filename = "%s_clusters.txt" % dirname
        
    print filename
    
    with open(filename, "w") as cluster_file:
        for frame in frames:
            cluster_file.write(frame + '\n')
            
