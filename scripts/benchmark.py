#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import datetime
from collections import defaultdict

results = """
01/31/2018 18:23:38	apinska	1744564 2ubq_bench_rigid	19:30:53	05:59:13	38	107020	1: Script or program ended with error
01/31/2018 19:58:05	apinska	1744565 2ubq_bench_flex	27:02:44	07:33:37	41	107219	1: Script or program ended with error

02/01/2018 06:48:49	apinska	1744600 2ubq_bench_rigid	19:55:24	06:03:59	38	107020	1: Script or program ended with error
02/01/2018 08:24:50	apinska	1744601 2ubq_bench_flex	27:35:54	07:39:54	39	107285	1: Script or program ended with error

02/01/2018 17:20:36	apinska	1744610 2ubq_bench_rigid	19:37:55	06:00:42	38	107020	1: Script or program ended with error
02/01/2018 18:59:32	apinska	1744611 2ubq_bench_flex	27:38:11	07:39:38	41	107285	134

02/02/2018 06:06:09	apinska	1744641 2ubq_bench_half	26:42:41	07:27:08	40	107218	1: Script or program ended with error
02/02/2018 06:42:41	apinska	1744642 2ubq_bench_all	29:11:36	08:03:37	40	107220	1: Script or program ended with error
02/02/2018 07:14:35	apinska	1744645 2ubq_bench_longtail	27:45:17	07:34:08	41	107283	1: Script or program ended with error

02/02/2018 15:58:16	apinska	1744647 2ubq_bench_half	26:35:50	07:26:17	40	107217	1: Script or program ended with error
02/02/2018 16:04:29	apinska	1744649 2ubq_bench_longtail	27:58:47	07:32:24	40	107217	1: Script or program ended with error
02/02/2018 16:35:21	apinska	1744648 2ubq_bench_all	29:33:28	08:03:18	41	107283	1: Script or program ended with error

02/03/2018 01:02:32	apinska	1744752 2ubq_bench_half	26:33:39	07:29:45	38	107219	134
02/03/2018 01:10:33	apinska	1744754 2ubq_bench_longtail	27:54:35	07:37:43	41	107285	1: Script or program ended with error
02/03/2018 01:35:56	apinska	1744753 2ubq_bench_all	29:15:51	08:03:07	40	107285	139: Segmentation fault / memory error

02/03/2018 16:10:12	apinska	1744766 2ubq_bench_rigid	19:47:11	06:02:42	38	107022	1: Script or program ended with error
02/03/2018 17:43:06	apinska	1744768 2ubq_bench_longtail	27:54:31	07:35:25	40	107219	1: Script or program ended with error
02/03/2018 17:46:17	apinska	1744767 2ubq_bench_flex	27:32:00	07:38:43	41	107285	1: Script or program ended with error

02/03/2018 23:49:49	apinska	1744769 2ubq_bench_half	27:32:23	07:39:36	40	107217	1: Script or program ended with error

02/04/2018 03:43:32	apinska	1744770 2ubq_bench_all	29:05:11	08:00:19	40	107217	1: Script or program ended with error
02/04/2018 08:01:04	apinska	1744771 2ubq_bench_rigid	19:38:05	06:00:23	36	107022	1: Script or program ended with error
02/04/2018 09:36:41	apinska	1744772 2ubq_bench_longtail	27:55:37	07:35:55	38	107219	1: Script or program ended with error

02/04/2018 19:12:16	apinska	1744773 2ubq_bench_half	26:22:31	07:21:55	40	107217	1: Script or program ended with error
02/04/2018 19:28:42	apinska	1744775 2ubq_bench_flex	28:01:33	07:37:28	41	107283	1: Script or program ended with error
02/04/2018 19:56:20	apinska	1744774 2ubq_bench_all	29:24:07	08:05:46	38	107220	1: Script or program ended with error
"""

def delta(s):
    h, m, s = (int(p) for p in s.split(':'))
    return datetime.timedelta(hours=h, minutes=m, seconds=s)

RESULT = re.compile('.*2ubq_bench_([a-z]+)	(\d+:\d+:\d+)	(\d+:\d+:\d+)	.*')

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

print("CPUTIME")

for k, v in cputimes.items():
    avg = sum(v, datetime.timedelta())/len(v)
    avg_hours = avg.days*24 + avg.seconds/(60*60)
    print("%s: %g hours (average of %d)" % (k, avg_hours, len(v)))
    averages[k]["cpu"] = avg_hours

print("WALLTIME")

for k, v in walltimes.items():
    avg = sum(v, datetime.timedelta())/len(v)
    avg_hours = avg.days*24 + avg.seconds/(60*60)
    print("%s: %g hours (average of %d)" % (k, avg_hours, len(v)))
    averages[k]["wall"] = avg_hours
    
for sim in ('rigid', 'flex', 'longtail', 'half', 'all'):
    print("%s: %.2f & %.2f" % (sim, averages[sim]["cpu"], averages[sim]["wall"]))
    
