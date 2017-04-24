#!/usr/bin/env python3
from caflib.Tools import geomlib
from caflib.Timing import timing
from caflib.Configure import before_templates, feature
from caflib.Tools.geomlib import Crystal
import vdwsets.vdwsets as vdw
import json
import os
from pathlib import Path


paths = [
    '<>/*/*/*/calc',
    '<>/*/*/{b3lyp*,lda,m06*,pbe*,scan*}',
    '*/*/*/<>/calc',
    '*/*/*/<b3lyp*,lda,m06*,pbe*,scan*>',
    '<>/*/*/mbd',
    '<>/*/*/dftd3',
    '<>/*/*/vv10/*/*/calc',
    '<>/*/*/vv10/*/*',
]
b_vv10_values = [4.5, 6.3, 8., 10., 13., 16., 19., 22.]

pseudo_root = Path(os.environ['ESPRESSO_PSEUDO'])


def taskgen(ctx, geom, dsname):
    with timing('taskgen'):
        master = ctx()
        for label, xc, basis in [
                ('lda', 'pw-lda', 'tight'),
                ('pbe', 'pbe', 'tight'),
                ('pbe0', 'pbe0', 'tight'),
                ('scan', 'dfauto scan', 'tight'),
                ('b3lyp', 'b3lyp', 'tight'),
                ('scan0', 'dfauto scan0', 'tight'),
                ('m06', 'm06', 'tight'),
                ('pbe-light', 'pbe', 'light'),
                ('scan-light', 'dfauto scan', 'light'),
                ('b3lyp-light', 'b3lyp', 'light'),
                ('pbe0-light', 'pbe0', 'light'),
                ('scan0-light', 'dfauto scan0', 'light'),
                ('m06-light', 'm06', 'light')
        ]:
            task = ctx(
                features='aims',
                templates=['geometry.in', (
                    'control.periodic.in' if hasattr(geom, 'lattice')
                    else 'control.in', 'control.in'
                )],
                geom=geom,
                charge=geom.metadata.get('charge') or 0,
                k_grid=geom.get_kgrid(0.8),
                output='hirshfeld' if label == 'pbe' else None,
                xc=xc,
                xc_pre=('m06-l', 10) if xc == 'm06' else None,
                basis=basis,
                aims='aims.master'
            ) + ctx.link('calc', ('run.out', 'aims.out')) + ctx(
                command='python3 process.py <aims.out >results.json',
                files='process.py'
            )
            if label == 'pbe':
                pbe = task
            task + ctx.link(label) + master
        inp = {
            'coords': geom.coords,
            'species': geom.species
        }
        if hasattr(geom, 'lattice'):
            inp['lattice'] = geom.lattice.tolist(),
            inp['k_grid'] = [int(x) for x in geom.get_kgrid(0.8)]
        pbe + ctx.link('scf', ('results.json', 'dft.json')) + ctx(
            command='run_mbd mbd.py',
            files='mbd.py',
            templates='input.json',
            input=json.dumps(inp, sort_keys=True)
        ) + ctx.link('mbd') + master
        if hasattr(geom, 'lattice'):
            tsk = ctx(
                command='time python3 dftd3-driver.py POSCAR -pbc -abc',
                files='dftd3-driver.py',
                texts={'POSCAR': geom.dumps('vasp')}
            )
        else:
            tsk = ctx(
                command='time python3 dftd3-driver.py geom.xyz -abc',
                files='dftd3-driver.py',
                texts={'geom.xyz': geom.dumps('xyz')}
            )
        tsk + ctx.link('dftd3') + master
        if dsname not in ['S22', 'X23', 'S66x8', 'S12L']:
            return master
        if hasattr(geom, 'lattice'):
            kwargs = {'k_grid': geom.get_kgrid(0.8), 'geom': geom}
        else:
            kwargs = {
                'k_grid': (1, 1, 1),
                'geom': Crystal.from_molecule(geom, padding=4.)
            }
        if 'charge' in geom.metadata:
            kwargs['other'] = f'\ntot_charge = {geom.metadata["charge"]}'
        vv10 = ctx()
        for b_vv10 in b_vv10_values:
            tsk = ctx()
            for xc_name, xc in [('vdw', 'rvv10'), ('base', 'sla+pw+rw86+pbc'), ('vdw2', 'rvv10')]:
                ctx(
                    features='qe',
                    templates=('qespresso.in', 'input.in'),
                    basis=30.,
                    pseudo='_ONCV_PBE-1.0.upf',
                    qe='qe.master' if xc_name == 'vdw2' else 'qe.8d90260',
                    b_vv=b_vv10,
                    xc=xc,
                    **kwargs
                ) + ctx.link('calc', ('run.out', 'qe.out')) + ctx(
                    command='python3 process_qe.py <qe.out >results.json',
                    files='process_qe.py'
                ) + ctx.link(xc_name) + tsk
            tsk + ctx.link(str(b_vv10)) + vv10
        vv10 + ctx.link('vv10') + master
    return master


def configure(ctx):
    ctx.load_tool('aims')
    geomlib.settings['eq_precision'] = 3
    for name, ds in vdw.get_all_datasets().items():
        ds.get_task(ctx, taskgen) * ctx.target(name)


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


@before_templates
@feature('qe')
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
        'k_grid': ' '.join(map(str, task.consume('k_grid'))),
    }
    attribs.update({
        k: v for k, v
        in ((attr, task.consume(attr)) for attr in 'basis xc b_vv C_vv'.split())
        if v
    })
    task.attrs.update(attribs)
    qe = task.consume('qe')
    task.attrs['command'] = f'QESPRESSO={qe} run_qe'
    task.inputs['geom.vasp'] = geom.dumps('vasp')
