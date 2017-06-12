import pandas as pd
from pathlib import Path
import numpy as np

root = '../calc/2017-02-13-graphene-flakes'
import sys
sys.path.append(root)
import caf  # noqa
from caflib.Cellar import Cellar

ev = 27.2107
kcal = 627.503


def collector():
    cellar = Cellar(f'{root}/.caf')
    tree = cellar.get_tree(objects=True)
    data = []
    for hashid, path in tree.glob('{dimers*,monomers}/*/*/*'):
        fragment, name, xc, basis = path.split('/')
        if 'strip' in name:
            continue
        if 'outputs' not in tree.objects[hashid]:
            continue
        filehash = tree.objects[hashid]['outputs']['run.out']
        with open(cellar.get_file(filehash)) as f:
            try:
                converged = False
                words = next(l for l in f if 'The structure contains' in l).split()
                natoms, nelec = int(words[3]), float(words[9])
                ncarbons = (int(nelec)-natoms)/5
                while True:
                    line = next(
                        l for l in f if
                        'Total energy uncorrected      ' in l or
                        'cycle converged' in l or
                        'HOMO-LUMO' in l
                    )
                    if 'converged' in line:
                        converged = True
                    elif 'HOMO' in line:
                        line = (w for w in line.split())
                        next(w for w in line if 'gap' in w)
                        gap = float(next(line))
                    else:
                        ene = float(line.split()[5])
            except StopIteration:
                pass
        if not converged:
            ene = np.NaN
        data.append((name, xc, basis, fragment, ene, nelec, ncarbons, gap))
    data = pd.DataFrame(
        data,
        columns='name xc basis fragment ene nelec ncarbons gap'.split()
    ).set_index('name xc basis fragment'.split(), verify_integrity=True).sort_index()
    return data


if __name__ == '__main__':
    data = collector()
    hdf_file = str(Path(__file__).parent/'../data/all-data.h5')
    with pd.HDFStore(hdf_file) as store:
        store['flakes'] = data
