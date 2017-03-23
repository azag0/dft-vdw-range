#!/usr/bin/env python3
import json
import sys

results = {}
f = sys.stdin
next(l for l in f if l == '  Self-consistency cycle converged.\n')
results['scf_energy'] = \
    float(next(l for l in f if 'Total energy uncorrected' in l).split()[5])
try:
    next(l for l in f if 'Performing Hirshfeld analysis' in l)
except StopIteration:
    pass
else:
    ratios = []
    for line in f:
        if not line.strip():
            break
        if 'Free atom volume' in line:
            free = float(line.split()[-1])
        elif 'Hirshfeld volume' in line:
            hirsh = float(line.split()[-1])
            ratios.append(hirsh/free)
    results['hirsh'] = ratios
json.dump(results, sys.stdout, sort_keys=True)
