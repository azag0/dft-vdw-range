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
    s22 = vdw.get_s22()
    data = []
    for key, cluster in s22.clusters.items():
        system, = key
        keypath = slugify('_'.join(map(str, key)))
        enes = {}
        params = {}
        for fragment, geomid in cluster.fragments.items():
            path = f'S22/{keypath}/{fragment}/pbe'
            if path not in tree:
                continue
            hashid = tree[path]
            results_file = 'results.json'
            filehash = tree.objects[hashid]['outputs'][results_file]
            with open(cellar.get_file(filehash)) as f:
                results = json.load(f)
            enes[fragment] = results['scf_energy']/ev*kcal
            params[fragment] = results['hirsh'], s22.geoms[geomid].species, s22.geoms[geomid].coords
        ene_int = cluster.get_int_ene(enes)
        data.append((system, ene_int, cluster.energies['ref'], params))
    return data


if __name__ == '__main__':
    data = collector()
    with open('s22.json', 'w') as f:
        json.dump(data, f)
    # hdf_file = str(Path(__file__).parent/'../data/s22.h5')
    # with pd.HDFStore(hdf_file) as store:
    #     for dfname, df in data.items():
    #         store[dfname] = df
