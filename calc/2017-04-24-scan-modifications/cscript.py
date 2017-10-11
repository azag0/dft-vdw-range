#!/usr/bin/env python3
from caflib.Tools import geomlib, aims
import vdwsets.vdwsets as vdw

aims._tags.append('xc_param')


def flatten(x):
    seqs = (tuple, list)
    if not isinstance(x, seqs):
        yield x
    else:
        for y in x:
            yield from flatten(y)


all_params = [
    [],
    ('scanx_k1', 0.03),
    ('scanx_k1', 0.1),
    ('scanx_k1', 0.4),
    ('scanx_c1x', 1.3),
    ('scanx_c1x', 2.),
    ('scanx_c2x', 0.3),
    ('scanx_c2x', 0.5),
    ('scanx_c2x', 1.4),
    ('scanx_dx', 2.0),
    ('scanx_dx', 1.5),
    ('scanx_dx', 1.),
    ('scanx_dx', 0.5),
    ('scanx_dx', 0.0),
    ('scanx_muak', 0.2195),
    [('scanx_k1', 0.4), ('scanx_muak', 0.2195)],
    [('scanx_hx0', 1.14), ('scanx_a1', 10.)],
    [('scanx_hx0', 1.2), ('scanx_a1', 3.38287)]
]


def taskgen(ctx, geom, dsname):
    master = ctx()
    for params in all_params:
        ctx(
            features='aims',
            templates='geometry.in control.in',
            geom=geom,
            charge=0,
            xc_param=params,
            xc='dfauto scan',
            basis='tight',
            aims_delink='aims.b4baa15' if geom.metadata.get('cp') else 'aims.1a57b68',
        ) + ctx.link('calc', ('run.out', 'aims.out')) + ctx(
            command='python3 process.py <aims.out >results.json',
            files='process.py'
        ) + ctx.link('_'.join(map(str, flatten(params))) if params else 'default') + master
    return master


def configure(ctx):
    ctx.load_tool('aims')
    geomlib.settings['eq_precision'] = 3
    for name, ds in vdw.get_all_datasets(include=['S22']).items():
        ds.get_task(ctx, taskgen) * ctx.target(name)
