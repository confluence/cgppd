#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import matplotlib.pyplot as plt

data = sys.stdin.read()
values = [float(f) for f in data.split('\n') if f]

plt.plot(values, 'bo')

if len(sys.argv) > 1:
    plt.savefig(sys.argv[1])
else:
    plt.show()
