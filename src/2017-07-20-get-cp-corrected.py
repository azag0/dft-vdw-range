import pandas as pd
import json

root = '../calc/2017-01-23-all-vdw-sets-3'
import sys
sys.path.append(root)
import caf  # noqa
from caflib.Cellar import Cellar
from caflib.Utils import slugify
import vdwsets.vdwsets as vdw

ev = 27.2107
kcal = 627.503

# ::>

cellar = Cellar(f'{root}/.caf')
tree = cellar.get_tree(objects=True)

# ::>

fragments = [
    ('scan', 'complex'),
    ('scan_-frag-0', 'fragment-1'),
    ('scan_-frag-1', 'fragment-2'),
]
data = []
for key, cluster in vdw.get_s66x8().clusters.items():
    keypath = slugify('_'.join(map(str, key)))
    system, dist = key
    ref = cluster.energies['ref']
    enes = {}
    for method, fragment in fragments:
        path = f'S66x8/{keypath}/complex/{method}'
        hashid = tree[path]
        if 'outputs' not in tree.objects[hashid]:
            continue
        results_file = 'results.json'
        filehash = tree.objects[hashid]['outputs'][results_file]
        with open(cellar.get_file(filehash)) as f:
            results = json.load(f)
        enes[fragment] = list(sorted(results)) if method in ['mbd', 'dftd3'] else results['scf_energy']
    if len(enes) != len(cluster.fragments):
        continue
    data.append((system, dist, cluster.get_int_ene(enes)/ev*kcal, ref))
pd.DataFrame(data, columns='system dist ene ref'.split())
