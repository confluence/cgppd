#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import os

vis_filename, output_desc, *cmd = sys.argv[1:]

NEW_MOLECULE = re.compile("mol new (.*\.pdb)")
CLUSTER = re.compile("mol drawframes top \d {(\d.*)}")
SIMNAME = re.compile("polyubiquitin(_ll)?/diubiquitin_(lys|met)_(\d+)_")

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

for dirname, frames_clusters in clusters.items():
    if os.path.isdir(dirname):
        filename = "%s/clusters_%s.txt" % (dirname, output_desc)
    else:
        filename = "%s_clusters_%s.txt" % (dirname, output_desc)
        
    total_frame_count = 0
    cluster_frame_counts = []
    representative_frames = []

    with open(filename, "w") as cluster_file:
        for frames_string in frames_clusters:
            frames = frames_string.split()
            
            total_frame_count += len(frames)
            cluster_frame_counts.append(len(frames))
            representative_frames.append(frames[0])
            
            if "save" in cmd:
                cluster_file.write(frames_string + '\n')
            

    if "save" in cmd:
        print("Saved to", filename)
    
    percentages = [int(100 * size / total_frame_count) for size in cluster_frame_counts]
    
    m = SIMNAME.search(dirname)
    simname = "diubq%s%s" % (m.group(3), "ll" if m.group(1) else "")
    
    if "perc" in cmd:
        print(simname, end=" ")
        for i, (size, perc) in enumerate(zip(cluster_frame_counts, percentages)):
            if perc >= 5:
                print("%d: %d/%d (%d%%)" % (i, size, total_frame_count, perc), end=" ")
        print()
                
    if "tcl" in cmd:
        print("# ", " ".join("%d%%" % perc for perc in percentages if perc >= 5))
        print("set cluster_frame(%s) { %s }" % (simname, " ".join([f for f, perc in zip(representative_frames, percentages) if perc >= 5])))
