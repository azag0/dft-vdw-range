#!/usr/bin/env python3
import json
import sys
import re


def get_seconds(time):
    time = re.split('([smh])', time)[:-1]
    return sum(float(num)*{'s': 1, 'm': 60, 'h': 3600}[unit] for num, unit in chunks(time, 2))


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


results = {}
f = sys.stdin
for line in f:
    if line[0] == '!':
        break
    if line.lstrip().startswith('Non-local correlation energy'):
        results['nlc'] = float(line.split()[-1])
results['energy'] = float(line.split()[4])
line = next(l for l in f if 'PWSCF        :' in l)
results['time'] = get_seconds(re.search(r'CPU +(\S.*) WALL', line).group(1))
json.dump(results, sys.stdout, sort_keys=True)
