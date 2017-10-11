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


def collector():
    cellar = Cellar(f'{root}/.caf')
    tree = cellar.get_tree(objects=True)
    s66 = vdw.get_s66x8()
    data = []
    for key, cluster in s66.clusters.items():
        system, dist = key
        keypath = slugify('_'.join(map(str, key)))
        enes = {}
        params = {}
        for fragment, geomid in cluster.fragments.items():
            path = f'S66x8/{keypath}/{fragment}/pbe'
            if path not in tree:
                continue
            hashid = tree[path]
            results_file = 'results.json'
            filehash = tree.objects[hashid]['outputs'][results_file]
            with open(cellar.get_file(filehash)) as f:
                results = json.load(f)
            enes[fragment] = results['scf_energy']/ev*kcal
            params[fragment] = results['hirsh'], s66.geoms[geomid].species, s66.geoms[geomid].coords
        ene_int = cluster.get_int_ene(enes)
        data.append((system, dist, ene_int, cluster.energies['ref'], params))
    return data


if __name__ == '__main__':
    data = collector()
    with open('s66.json', 'w') as f:
        json.dump(data, f)
