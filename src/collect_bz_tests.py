import pandas as pd
from pathlib import Path
import re

root = '../calc/2017-03-24-benzene-dimer-qe-tests/'
import sys
sys.path.append(root)
import caf  # noqa
from caflib.Cellar import Cellar

ev = 27.2107
kcal = 627.503


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def get_seconds(time):
    time = re.split('([smh])', time)[:-1]
    return sum(float(num)*{'s': 1, 'm': 60, 'h': 3600}[unit] for num, unit in chunks(time, 2))


def collector():
    cellar = Cellar(f'{root}/.caf')
    tree = cellar.get_tree(objects=True)
    data = []
    for hashid, path in tree.glob('**'):
        if 'outputs' not in tree.objects[hashid]:
            continue
        keys = path.split('/')
        filehash = tree.objects[hashid]['outputs']['run.out']
        with open(cellar.get_file(filehash)) as f:
            line = next(l for l in f if l[0] == '!')
            ene = float(line.split()[4])*kcal/2  # from Ry
            line = next(l for l in f if 'PWSCF        :' in l)
            time = get_seconds(re.search(r'CPU +(\S.*) WALL', line).group(1))
        data.append((*keys, ene, time))
    data = pd.DataFrame(
        data,
        columns='system padding basis pseudo xc subsystem ene time'.split()
    ).set_index('system padding basis pseudo xc subsystem'.split(), verify_integrity=True) \
        .sort_index()
    return data


if __name__ == '__main__':
    data = collector()
    hdf_file = str(Path(__file__).parent/'../data/all-data.h5')
    with pd.HDFStore(hdf_file) as store:
        store['bz_tests'] = data
