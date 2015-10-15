#!/usr/bin/env python

# Plot histogram produced by g_cluster

import sys
import re
import matplotlib.pyplot as plt

try:
    filename = sys.argv[1]
except IOError:
    sys.exit("No file provided.")

headers = []
items = []

with open(filename, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        elif line.startswith("@"):
            headers.append(line)
            continue
            
        x, y = line.split()
        items.append((float(x), int(y)))
        
header_str = "".join(headers)
title, xlabel, ylabel = re.search('title "(.*?)".*xaxis  label "(.*?)".*yaxis  label "(.*?)"', header_str, re.MULTILINE|re.DOTALL).groups()

simdir = re.search('(diubiquitin[^/]+)', filename).group(1)
title = "%s - %s" % (title, simdir)
            
x, y = zip(*items)

plt.bar(x, y, width=0.01)
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xlim([0, 2.5])
plt.ylim([0, 1400000])

plt.show()
