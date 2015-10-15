#!/usr/bin/env python

# Stupid demo for investigating maths behind distance histograms

import sys
import matplotlib.pyplot as plt
from collections import Counter

things = [float(arg) for arg in sys.argv[1:]]
pairs = [(t1, t2) for i, t1 in enumerate(things) for t2 in things[i+1:]]
dists = [abs(t1 - t2) for (t1, t2) in pairs]
rounded_dists = [int(d) for d in dists]
hist = Counter(rounded_dists)

x, y = zip(*sorted(hist.items()))
plt.bar(x, y)
plt.show()
