#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import os

vis_filename, output_desc = sys.argv[1:3]

NEW_MOLECULE = re.compile("mol new (.*\.pdb)")
CLUSTER = re.compile("mol drawframes top \d {(\d.*)}")

clusters = {}
current_cluster = None

with open(vis_filename, "r") as vis_file:
    for line in vis_file:
        nm = NEW_MOLECULE.match(line)
        if nm:
            current_cluster = os.path.dirname(nm.group(1))
            clusters[current_cluster] = []
            continue
            
        cl = CLUSTER.match(line)
        if cl:
            clusters[current_cluster].append(cl.group(1))

for dirname, frames in clusters.iteritems():
    print dirname
    if os.path.isdir(dirname):
        filename = "%s/clusters_%s.txt" % (dirname, output_desc)
    else:
        filename = "%s_clusters_%s.txt" % (dirname, output_desc)
        
    print filename
    
    with open(filename, "w") as cluster_file:
        for frame in frames:
            cluster_file.write(frame + '\n')
            
