import pandas as pd
import gc
import json
from collections import defaultdict
from itertools import product
from pathlib import Path

root = '../calc/2017-01-23-all-vdw-sets-3'
import sys
sys.path.append(root)
import caf  # noqa
from caflib.Cellar import Cellar
from caflib.Utils import slugify
import vdwsets.vdwsets as vdw
from cscript import b_vv10_values

ev = 27.2107
kcal = 627.503


def get_some(d, *keys):
    for key in keys:
        if key in d:
            return d[key]


def collector(methods=None, **kwargs):
    cellar = Cellar(f'{root}/.caf')
    tree = cellar.get_tree(objects=True)
    methods = methods or [p.split('/')[-1] for p in tree.dglob('S22/*/*/<>')]
    data = defaultdict(list)
    for name, ds in vdw.get_all_datasets(**kwargs).items():
        for key, cluster in ds.clusters.items():
            def get_int_ene(enes, cp=False):
                try:
                    ene_int = cluster.get_int_ene(enes)
                except KeyError:
                    return
                if cp and isinstance(ene_int, float):
                    ene_int = (ene_int, None)
                return ene_int

            fragments_nocp = [f for f in cluster.fragments if not f.endswith('-cp') and not f.endswith('-nocp')]
            keypath = slugify('_'.join(map(str, key)))
            system, dist = key if len(key) == 2 else (key[0], 1.)
            refs = cluster.energies
            ref = refs['ref'] if 'ref' in refs else refs['ref-grimme']
            for method in {*methods} - {'vv10'}:
                enes = {}
                for fragment in cluster.fragments:
                    path = f'{name}/{keypath}/{fragment}/{method}'
                    if path not in tree:
                        continue
                    hashid = tree[path]
                    if 'outputs' not in tree.objects[hashid]:
                        continue
                    results_file = 'mbd.json' if method == 'mbd' else 'results.json'
                    filehash = tree.objects[hashid]['outputs'][results_file]
                    with open(cellar.get_file(filehash)) as f:
                        results = json.load(f)
                    enes[fragment] = list(sorted(results)) if method in ['mbd', 'dftd3'] else results['scf_energy']
                if get_int_ene(dict.fromkeys(enes, 0)) is None:
                    continue
                if method == 'mbd':
                    for i, (beta, a, _) in enumerate(list(enes.values())[0]):
                        data['mbd'].append((
                            name, system, dist, a, beta, cluster.get_int_ene(
                                {f: enes[f][i][2] for f in fragments_nocp}
                            )*kcal,
                        ))
                    data['mbd'].append((name, system, dist, 6., 10., 0.))
                elif method == 'dftd3':
                    for i, (damping, *params, _, _) in enumerate(list(enes.values())[0]):
                        (data['dftd3_0'] if damping == 'zero' else data['dftd3_bj']).append((
                            name, system, dist, *params,
                            get_int_ene(
                                {f: enes[f][i][-2] for f in fragments_nocp}
                            )*kcal,
                            get_int_ene(
                                {f: enes[f][i][-1] for f in fragments_nocp}
                            )*kcal
                        ))
                    data['dftd3_bj'].append((name, system, dist, 1., 50., 1., 0., 0.))
                    data['dftd3_0'].append((name, system, dist, 10., 10., 1., 0., 0.))
                else:
                    ene_int = get_int_ene(enes, cp=True)
                    for ene, cp in zip(ene_int, [False, True]):
                        if ene is None:
                            continue
                        data['scf'].append((
                            name, system, dist, method, cp,
                            ene/ev*kcal, ref,
                            len(ds.geoms[get_some(cluster.fragments, 'complex', 'dimer', 'crystal', 'ABC')])
                        ))
            if name not in ['S22', 'X23', 'S66x8', 'S12L']:
                continue
            enes = defaultdict(dict)
            for fragment in fragments_nocp:
                for xc, b_vv10 in product(['base', 'vdw', 'vdw2'], b_vv10_values):
                    path = f'{name}/{keypath}/{fragment}/vv10/{b_vv10}/{xc}'
                    hashid = tree[path]
                    if 'outputs' not in tree.objects[hashid]:
                        continue
                    filehash = tree.objects[hashid]['outputs']['results.json']
                    with open(cellar.get_file(filehash)) as f:
                        results = json.load(f)
                    enes[b_vv10, xc][fragment] = results['energy']
                    data['vv10_raw'].append((
                        name, system, dist, b_vv10, xc, fragment, results['energy']*kcal/2
                    ))
                    if xc == 'vdw2':
                        enes[b_vv10, 'nlc'][fragment] = results['nlc']
                        data['vv10_raw'].append((
                            name, system, dist, b_vv10, 'nlc', fragment, results['nlc']*kcal/2
                        ))
            for (b_vv10, xc), enes in enes.items():
                data['vv10'].append((
                    name, system, dist, b_vv10, xc, cluster.get_int_ene(enes)*kcal/2
                ))
    columns = {
        'scf': 'name system dist xc cp|ene ref natoms',
        'mbd': 'name system dist a beta|ene',
        'dftd3_0': 'name system dist sr6 sr8 s8|ene ene3',
        'dftd3_bj': 'name system dist a1 a2 s8|ene ene3',
        'vv10': 'name system dist b xc|ene',
        'vv10_raw': 'name system dist b xc fragment|ene',
    }
    for dfname, rows in data.items():
        idx_cols, val_cols = (cols.split() for cols in columns[dfname].split('|'))
        data[dfname] = pd.DataFrame(rows, columns=idx_cols + val_cols) \
            .set_index(idx_cols, verify_integrity=True).sort_index()
    gc.collect()
    return data


if __name__ == '__main__':
    data = collector()
    hdf_file = str(Path(__file__).parent/'../data/all-data.h5')
    with pd.HDFStore(hdf_file) as store:
        for dfname, df in data.items():
            store[dfname] = df
