#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import datetime
from collections import defaultdict

benchmark_runs = {
    "on GPU003; 10 CPUs each, 4 GPUs each" : """
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
""",
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
        
