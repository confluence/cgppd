#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import datetime
from collections import defaultdict

benchmark_runs = {
    #"on GPU004; 20 CPUs each, 2 GPUs each" : """
#02/08/2018 13:41:56	apinska	1745269 2ubq_bench_all	75:55:15	04:06:47	179	168703	0: OK
#02/08/2018 18:45:43	apinska	1745270 2ubq_bench_flex	95:10:53	05:00:33	183	168965	0: OK
#02/08/2018 23:37:58	apinska	1745271 2ubq_bench_half	90:23:47	04:47:44	177	168571	0: OK
#02/09/2018 04:36:14	apinska	1745272 2ubq_bench_longtail	93:39:18	04:56:02	180	168833	0: OK
#02/09/2018 11:26:15	apinska	1745273 2ubq_bench_rigid	132:13:38	06:46:03	175	168178	0: OK

#02/09/2018 19:07:47	apinska	1745489 2ubq_bench_all	76:10:36	04:07:32	179	168703	0: OK
#02/10/2018 00:10:45	apinska	1745490 2ubq_bench_flex	95:15:08	05:00:31	182	168964	0: OK
#02/10/2018 04:59:48	apinska	1745491 2ubq_bench_half	89:39:11	04:44:40	177	168571	0: OK
#02/10/2018 09:54:26	apinska	1745492 2ubq_bench_longtail	92:53:16	04:54:14	180	168833	0: OK
#02/10/2018 16:42:52	apinska	1745493 2ubq_bench_rigid	132:38:41	06:47:45	173	168178	0: OK
#""",
    "on GPU003; 10 CPUs each, 4 GPUs each" : """
02/10/2018 21:16:03	apinska	1745527 2ubq_bench_all	21:18:46	02:25:52	332	194251	0: OK
02/10/2018 23:50:27	apinska	1745528 2ubq_bench_flex	21:53:25	02:30:18	331	194250	0: OK
02/11/2018 02:11:15	apinska	1745529 2ubq_bench_half	18:12:45	02:16:08	331	194251	0: OK
02/11/2018 04:36:03	apinska	1745530 2ubq_bench_longtail	20:02:59	02:20:54	331	194250	0: OK
02/11/2018 07:17:04	apinska	1745531 2ubq_bench_rigid	18:13:07	02:36:55	328	194250	0: OK
""",
    #"on GPU003; 12 CPUs each, 4 GPUs each" : """
#""",
}

def delta(s):
    h, m, s = (int(p) for p in s.split(':'))
    return datetime.timedelta(hours=h, minutes=m, seconds=s)

RESULT = re.compile('.*2ubq_bench_([a-z]+)	(\d+:\d+:\d+)	(\d+:\d+:\d+)	.*')

for params, results in benchmark_runs.items():
    
    print()
    print(params)

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
        averages[k]["cpu"] = avg_hours

    print("\nWALLTIME")

    for k, v in walltimes.items():
        avg = sum(v, datetime.timedelta())/len(v)
        avg_hours = avg.days*24 + avg.seconds/(60*60)
        print("%s: %g hours (average of %d)" % (k, avg_hours, len(v)))
        averages[k]["wall"] = avg_hours

    print("\nLATEX")
        
    for sim in ('rigid', 'flex', 'longtail', 'half', 'all'):
        print("%s: %.2f & %.2f" % (sim, averages[sim]["cpu"], averages[sim]["wall"]))
        
