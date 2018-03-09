#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import datetime
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

results = """
02/10/2018 21:16:03	apinska	1745527 2ubq_bench_all	21:18:46	02:25:52	332	194251	0: OK
02/10/2018 23:50:27	apinska	1745528 2ubq_bench_flex	21:53:25	02:30:18	331	194250	0: OK
02/11/2018 02:11:15	apinska	1745529 2ubq_bench_half	18:12:45	02:16:08	331	194251	0: OK
02/11/2018 04:36:03	apinska	1745530 2ubq_bench_longtail	20:02:59	02:20:54	331	194250	0: OK
02/11/2018 07:17:04	apinska	1745531 2ubq_bench_rigid	18:13:07	02:36:55	328	194250	0: OK

02/11/2018 15:00:27	apinska	1745566 2ubq_bench_all	21:15:05	02:25:20	332	194251	0: OK
02/11/2018 17:35:10	apinska	1745567 2ubq_bench_flex	21:49:16	02:30:03	331	194250	0: OK
02/11/2018 19:57:09	apinska	1745568 2ubq_bench_half	18:27:39	02:17:00	331	194251	0: OK
02/11/2018 22:20:54	apinska	1745569 2ubq_bench_longtail	20:03:26	02:20:29	331	194250	0: OK
02/12/2018 00:59:03	apinska	1745570 2ubq_bench_rigid	18:19:29	02:33:46	328	194250	0: OK
"""

def delta(s):
    h, m, s = (int(p) for p in s.split(':'))
    return datetime.timedelta(hours=h, minutes=m, seconds=s)

RESULT = re.compile('.*2ubq_bench_([a-z]+)	(\d+:\d+:\d+)	(\d+:\d+:\d+)	.*')
SIMS = ('rigid', 'flex', 'longtail', 'half', 'all')

cputimes = defaultdict(list)
walltimes = defaultdict(list)

for result in results.split('\n'):
    if not result.strip():
        continue
    
    m = RESULT.search(result)
    name, cputime, walltime = m.groups()
    
    cputimes[name].append(delta(cputime))
    walltimes[name].append(delta(walltime))
    
averages = defaultdict(dict)

print("\nCPUTIME")

for k, v in cputimes.items():
    avg = sum(v, datetime.timedelta())/len(v)
    avg_hours = avg.days*24 + avg.seconds/(60*60)
    print("%s: %g hours (average of %d)" % (k, avg_hours, len(v)))
    averages[k]["cpu"] = avg_hours / 10

print("\nWALLTIME")

for k, v in walltimes.items():
    avg = sum(v, datetime.timedelta())/len(v)
    avg_hours = avg.days*24 + avg.seconds/(60*60)
    print("%s: %g hours (average of %d)" % (k, avg_hours, len(v)))
    averages[k]["wall"] = avg_hours

print("TIMERS")

timers = defaultdict(dict)

for simname in SIMS:
    copy = []
    bonded = []
    nonbonded = []

    for filename in glob.glob("output/benchmark_logs/gpu3_10_4/*%s*" % simname):
        
        with open(filename) as logfile:
            text = logfile.read()
        
        copy.append(sum(float(v) for v in re.findall('(\d+\.\d+) \d+\.\d+ (?:Replica to GPU|Update replica on GPU|Update molecule on GPU)', text)))
        
        bonded.append(sum(float(v) for v in re.findall('(\d+\.\d+) \d+\.\d+ Bonded potential', text)))
        
        nonbonded.append(sum(float(v) for v in re.findall('(\d+\.\d+) \d+\.\d+ Kernel', text)))
        
    timers["copy"][simname] = sum(copy)/len(copy)
    timers["bonded"][simname] = sum(bonded)/len(bonded)
    timers["nonbonded"][simname] = sum(nonbonded)/len(nonbonded)
    
print(timers)

plt.figure()
#plt.xticks([0, 0.0526, 0.1711, 0.5, 1], ('rigid', 'short', 'long', 'half', 'all'))
#plt.plot(vals("cpu"), 'ro')
plt.plot([0, 0.0526, 0.1711, 0.5, 1], [averages[sim]["wall"] for sim in SIMS], 'bo-')
plt.xlabel("Fraction of residues in flexible linker")
plt.ylabel("Simulation running time (hours)")
fig = plt.gcf()
fig.savefig("walltime.svg", format='svg')

plt.figure()
plt.plot([0, 0.0526, 0.1711, 0.5, 1], [timers["copy"][sim] for sim in SIMS], 'mo-', label="data transfer")
plt.plot([0, 0.0526, 0.1711, 0.5, 1], [timers["bonded"][sim] for sim in SIMS], 'rs-', label="bonded potential")
plt.plot([0, 0.0526, 0.1711, 0.5, 1], [timers["nonbonded"][sim] for sim in SIMS], 'g^-', label="non-bonded potential")
plt.xlabel("Fraction of residues in flexible linker")
plt.ylabel("Total execution time (ms)")
plt.legend()
fig = plt.gcf()
fig.savefig("timers.svg", format='svg')
