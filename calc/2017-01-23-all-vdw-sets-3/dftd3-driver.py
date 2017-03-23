#!/usr/bin/env python3
from itertools import product
import subprocess as sp
import json
import os
import shutil
import sys
import numpy as np
import time
from tempfile import TemporaryDirectory
from multiprocessing import Process, Queue
import re


kcal = 627.503

sr6s = [0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.7, 2., 2.5, 3.]
sr8s = [0.6, 0.8, 1., 1.2, 1.4]
s8s = [-.3, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.7, 2., 2.5, 3.]
a1s = np.linspace(0.2, 1, 17)
a2s = np.linspace(0, 10, 31)

paramfile = '.dftd3par.local'


def worker(queuein, queueout):
    timer = 0.
    with TemporaryDirectory() as tmpdir:
        shutil.copy(sys.argv[1], tmpdir)
        while True:
            task = queuein.get()
            if task is None:
                break
            damping, *params = task
            assert damping in ['zero', 'BJ']
            with open(os.path.join(tmpdir, paramfile), 'w') as f:
                if damping == 'zero':
                    sr6, sr8, s8 = params
                    f.write(f'1. {sr6} {s8} {sr8} 14 3')
                elif damping == 'BJ':
                    a1, a2, s8 = params
                    f.write(f'1. {a1} {s8} {a2} 0 4')
            timer -= time.time()
            proc = sp.run(['dftd3', *sys.argv[1:]], stdout=sp.PIPE, cwd=tmpdir)
            timer += time.time()
            for line in proc.stdout.decode().split('\n'):
                line = line.lstrip()
                if line.startswith('Edisp'):
                    m = re.search(r'Edisp /kcal,au(?:,eV)?: *-?\d+\.\d+ *([- ]\d+\.\d+)', line)
                    if m:
                        ene = float(m.group(1))
                    else:
                        print(line)
                        ene = np.nan
                elif line.startswith('E6(ABC)'):
                    ene3 = float(line.split()[-1])/kcal
                    break
            queueout.put((*task, ene, ene3))
    queueout.put(timer)


queuein = Queue()
queueout = Queue()
nworkers = os.cpu_count()
pool = [
    Process(target=worker, args=(queuein, queueout)) for _ in range(nworkers)
]
for process in pool:
    process.start()
ntask = 0
for sr6, sr8, s8 in product(sr6s, sr8s, s8s):
    queuein.put(('zero', sr6, sr8, s8))
    ntask += 1
for a1, a2, s8 in product(a1s, a2s, s8s):
    queuein.put(('BJ', a1, a2, s8))
    ntask += 1
for _ in pool:
    queuein.put(None)
enes = []
timer = 0
while ntask > -nworkers:
    result = queueout.get()
    if isinstance(result, float):
        timer += result
    else:
        enes.append(result)
    ntask -= 1
for process in pool:
    process.join()
with open('results.json', 'w') as f:
    json.dump(enes, f)
print(timer)
