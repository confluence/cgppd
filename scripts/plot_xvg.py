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
            
x, y = zip(*items)

print x
print y

plt.bar(x, y, width=0.01)
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.show()
