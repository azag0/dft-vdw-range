#!/usr/bin/env python3
from caflib.Tools import geomlib
from caflib.Timing import timing
import vdwsets.vdwsets as vdw
import json


paths = [
    '<>/*/*/*/calc',
    '<>/*/*/{b3lyp*,lda,m06*,pbe*,scan*}',
    '*/*/*/<>/calc',
    '*/*/*/<b3lyp*,lda,m06*,pbe*,scan*>',
    '<>/*',
    '<>/*/*/mbd',
    '<>/*/*/dftd3',
]


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
    return master


def configure(ctx):
    ctx.load_tool('aims')
    geomlib.settings['eq_precision'] = 3
    for name, ds in vdw.get_all_datasets().items():
        ds.get_task(ctx, taskgen) * ctx.target(name)
