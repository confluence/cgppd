#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import datetime
from collections import defaultdict

benchmark_runs = {
    #"on GPU002; 4 CPUs each, 1 GPU each" : """
#01/31/2018 18:23:38	apinska	1744564 2ubq_bench_rigid	19:30:53	05:59:13	38	107020	1: Script or program ended with error
#01/31/2018 19:58:05	apinska	1744565 2ubq_bench_flex	27:02:44	07:33:37	41	107219	1: Script or program ended with error

#02/01/2018 06:48:49	apinska	1744600 2ubq_bench_rigid	19:55:24	06:03:59	38	107020	1: Script or program ended with error
#02/01/2018 08:24:50	apinska	1744601 2ubq_bench_flex	27:35:54	07:39:54	39	107285	1: Script or program ended with error

#02/01/2018 17:20:36	apinska	1744610 2ubq_bench_rigid	19:37:55	06:00:42	38	107020	1: Script or program ended with error
#02/01/2018 18:59:32	apinska	1744611 2ubq_bench_flex	27:38:11	07:39:38	41	107285	134

#02/02/2018 06:06:09	apinska	1744641 2ubq_bench_half	26:42:41	07:27:08	40	107218	1: Script or program ended with error
#02/02/2018 06:42:41	apinska	1744642 2ubq_bench_all	29:11:36	08:03:37	40	107220	1: Script or program ended with error
#02/02/2018 07:14:35	apinska	1744645 2ubq_bench_longtail	27:45:17	07:34:08	41	107283	1: Script or program ended with error

#02/02/2018 15:58:16	apinska	1744647 2ubq_bench_half	26:35:50	07:26:17	40	107217	1: Script or program ended with error
#02/02/2018 16:04:29	apinska	1744649 2ubq_bench_longtail	27:58:47	07:32:24	40	107217	1: Script or program ended with error
#02/02/2018 16:35:21	apinska	1744648 2ubq_bench_all	29:33:28	08:03:18	41	107283	1: Script or program ended with error

#02/03/2018 01:02:32	apinska	1744752 2ubq_bench_half	26:33:39	07:29:45	38	107219	134
#02/03/2018 01:10:33	apinska	1744754 2ubq_bench_longtail	27:54:35	07:37:43	41	107285	1: Script or program ended with error
#02/03/2018 01:35:56	apinska	1744753 2ubq_bench_all	29:15:51	08:03:07	40	107285	139: Segmentation fault / memory error

#02/03/2018 16:10:12	apinska	1744766 2ubq_bench_rigid	19:47:11	06:02:42	38	107022	1: Script or program ended with error
#02/03/2018 17:43:06	apinska	1744768 2ubq_bench_longtail	27:54:31	07:35:25	40	107219	1: Script or program ended with error
#02/03/2018 17:46:17	apinska	1744767 2ubq_bench_flex	27:32:00	07:38:43	41	107285	1: Script or program ended with error

#02/03/2018 23:49:49	apinska	1744769 2ubq_bench_half	27:32:23	07:39:36	40	107217	1: Script or program ended with error

#02/04/2018 03:43:32	apinska	1744770 2ubq_bench_all	29:05:11	08:00:19	40	107217	1: Script or program ended with error
#02/04/2018 08:01:04	apinska	1744771 2ubq_bench_rigid	19:38:05	06:00:23	36	107022	1: Script or program ended with error
#02/04/2018 09:36:41	apinska	1744772 2ubq_bench_longtail	27:55:37	07:35:55	38	107219	1: Script or program ended with error

#02/04/2018 19:12:16	apinska	1744773 2ubq_bench_half	26:22:31	07:21:55	40	107217	1: Script or program ended with error
#02/04/2018 19:28:42	apinska	1744775 2ubq_bench_flex	28:01:33	07:37:28	41	107283	1: Script or program ended with error
#02/04/2018 19:56:20	apinska	1744774 2ubq_bench_all	29:24:07	08:05:46	38	107220	1: Script or program ended with error
#""",
    #"on GPU004; 4 CPUs each, sharing 2 GPUs" : """
#02/08/2018 04:26:10	apinska	1745054 2ubq_bench_flex	68:35:50	25:15:55	174	168243	0: OK
#02/08/2018 05:44:18	apinska	1745055 2ubq_bench_longtail	72:40:48	26:34:02	174	168243	0: OK
#02/08/2018 06:24:21	apinska	1745053 2ubq_bench_rigid	69:51:38	27:14:07	171	168111	0: OK
#02/08/2018 06:34:07	apinska	1745057 2ubq_bench_all	73:17:30	27:23:49	174	168243	0: OK
#02/08/2018 06:34:54	apinska	1745056 2ubq_bench_half	71:17:35	27:24:37	174	168243	0: OK
#""",
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
        
