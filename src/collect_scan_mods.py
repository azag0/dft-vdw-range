import pandas as pd
from pathlib import Path
import json
from math import nan

root = '../calc/2017-04-24-scan-modifications'
import sys
sys.path.append(root)
import caf  # noqa
from caflib.Cellar import Cellar
from caflib.Utils import slugify
import vdwsets.vdwsets as vdw
from cscript import flatten, all_params

ev = 27.2107
kcal = 627.503


def collector():
    cellar = Cellar(f'{root}/.caf')
    tree = cellar.get_tree(objects=True)
    data = []
    for name, ds in vdw.get_all_datasets(include=['S22']).items():
        for key, cluster in ds.clusters.items():
            keypath = slugify('_'.join(map(str, key)))
            system = key[0]
            ref = cluster.energies['ref']
            for params in all_params:
                param_name = '_'.join(map(str, flatten(params))) if params else 'default'
                enes = {}
                for fragment in cluster.fragments:
                    path = f'{name}/{keypath}/{fragment}/{param_name}'
                    hashid = tree[path]
                    if 'outputs' not in tree.objects[hashid]:
                        enes[fragment] = nan
                        continue
                    filehash = tree.objects[hashid]['outputs']['results.json']
                    with open(cellar.get_file(filehash)) as f:
                        results = json.load(f)
                    enes[fragment] = results['scf_energy']
                data.append((
                    name, system, ' '.join(map(str, flatten(params))),
                    cluster.get_int_ene(enes)[1]/ev*kcal, ref
                ))
    data = pd.DataFrame(
        data,
        columns='name system params ene ref'.split()
    ).set_index('name system params'.split(), verify_integrity=True).sort_index()
    return data


if __name__ == '__main__':
    data = collector()
    hdf_file = str(Path(__file__).parent/'../data/all-data.h5')
    with pd.HDFStore(hdf_file) as store:
        store['scan_mods'] = data
