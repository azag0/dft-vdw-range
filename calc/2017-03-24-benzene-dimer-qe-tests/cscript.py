import shutil
import os
from pathlib import Path
from itertools import product

from caflib.Tools import geomlib
from caflib.Configure import before_templates


pseudo_root = Path(os.environ['ESPRESSO_PSEUDO'])


def configure(ctx):
    crystal = geomlib.readfile('crystal.aims')
    molecule = geomlib.readfile('molecule.aims')
    other = []
    for geomfile in ['bz-dimer-long.aims', 'supra.aims']:
        cluster = geomlib.readfile(geomfile)
        frags = cluster.get_fragments()
        system = Path(geomfile).stem
        other.extend([
            (system, 'complex', cluster),
            (system, 'frag-1a', frags[0]),
            (system, 'frag-2a', frags[1]),
        ])
    for padding, (xc_name, xc), basis, (pp_name, pp) in product(
            [2., 3., 4., 5., 6., 7.],
            [('vv10', 'rvv10'), ('base', 'sla+pw+rw86+pbc')],
            [20., 30., 50., 70.],
            [
                ('NC', '_ONCV_PBE-1.0.upf'),
                ('NC-SL', '.pbe-mt_fhi.UPF'),
                ('PAW', '.pbe-*kjpaw?psl*.UPF'),
                ('US', '.pbe-*rrkjus?psl*.UPF'),
                ('US-LDA', '.pw91*van*.UPF'),
            ],
    ):
        if padding > 5 and not (basis == 30. and pp_name == 'NC'):
            continue
        geom = geomlib.readfile('bz-dimer.xyz')
        geom_pbc = geomlib.Crystal.from_molecule(geom, padding)
        frags = geom_pbc.get_fragments()
        molecule_pbc = geomlib.Crystal.from_molecule(molecule, padding)
        all_geoms = [
            ('benzene', 'complex', geom_pbc),
            ('benzene', 'frag-1', geomlib.Crystal(geom_pbc.lattice, frags[0].atoms)),
            ('benzene', 'frag-2', geomlib.Crystal(geom_pbc.lattice, frags[1].atoms)),
            ('benzene', 'frag-1a', geomlib.Crystal.from_molecule(frags[0], padding)),
            ('benzene', 'frag-2a', geomlib.Crystal.from_molecule(frags[1], padding)),
            ('benzene', 'crystal', crystal),
            ('benzene', 'molecule', molecule_pbc),
        ]
        if basis == 30. and pp_name == 'NC':
            all_geoms.extend(
                (system, name, geomlib.Crystal.from_molecule(geom, padding))
                for system, name, geom in other
            )
        for system, name, geom in all_geoms:
            ctx(
                features=qespresso,
                templates=('qespresso.in', 'input.in'),
                basis=basis,
                pseudo=pp,
                k_grid=geom.get_kgrid(0.8) if name == 'crystal' else (1, 1, 1),
                qe_delink='qe.8d90260',
                geom=geom,
                xc=xc,
            ) * ctx.target(f'{system}/{padding}/{basis}/{pp_name}/{xc_name}/{name}')


class Cache:
    qe = {}
    pseudo = {}


def get_pp(elem, pseudo):
    if (elem, pseudo) in Cache.pseudo:
        return Cache.pseudo[elem, pseudo]
    candidates = [p.name for p in pseudo_root.glob(elem + pseudo)]
    if len(candidates) > 1:
        candidates = [c for c in candidates if 'ak' not in c]
    if len(candidates) != 1:
        raise RuntimeError((elem, pseudo, candidates))
    Cache.pseudo[elem, pseudo] = candidates[0]
    print(f'{elem}{pseudo} is {candidates[0]}')
    return candidates[0]


def get_qe(qe):
    if qe in Cache.qe:
        return Cache.qe[qe]
    exe = Path(shutil.which(qe))
    if exe.is_symlink():
        exe = os.readlink(exe)
    else:
        exe = qe
    Cache.qe[qe] = exe
    print(f'QE {qe} is {exe}')
    return exe


@before_templates
def qespresso(task):
    geom = task.consume('geom')
    pseudo = task.consume('pseudo')
    attribs = {
        'nat': len(geom),
        'ntyp': len(set(a.symbol for a in geom)),
        'species': '\n'.join(sorted(set(
            f'{a.symbol} {a.mass} {get_pp(a.symbol, pseudo)}' for a in geom
        ))),
        'coords': '\n'.join(' '.join(map(str, [a.symbol, *a.xyz])) for a in geom),
        'lattice': '\n'.join(' '.join(map(str, vec)) for vec in geom.lattice),
        'k_grid': ' '.join(str(int(x)) for x in task.consume('k_grid')),
    }
    attribs.update({
        k: v for k, v
        in ((attr, task.consume(attr)) for attr in 'basis xc b_vv C_vv'.split())
        if v
    })
    task.attrs.update(attribs)
    qe = get_qe(task.consume('qe_delink'))
    task.attrs['command'] = f'QESPRESSO={qe} run_qe'
    task.inputs['geom.vasp'] = geom.dumps('vasp')
