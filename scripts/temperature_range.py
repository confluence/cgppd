#!/usr/bin/env python

import sys

t_min = float(sys.argv[1])
t_max = float(sys.argv[2])
replicas = int(sys.argv[3])

t_geom = (t_max / t_min) ** (1.0 / (replicas - 1))

temps = [t_min * (t_geom ** i) for i in range(replicas)]

print "%.1f - %.1f / %d: %s" % (t_min, t_max, replicas, " ".join("%.1f" % t for t in temps))
