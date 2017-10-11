from caflib.Tools.geomlib import Atom, Molecule, Crystal
import numpy as np
from numpy import sqrt, cos, sin, pi
from numpy.linalg import norm
import re

from caflib.Configure import feature

paths = ['*/*/<>/<>', 'relax/*']

a = 2.461
dCH = 1.1
CC = 3.35


def configure(ctx):
    info = {
        'graphene': {
            'k_grid': (16, 16, 1),
            'vacuum': (False, False, True)
        }
    }
    objects = {
        'benzene': [(0, 0)],
        'naphtalene': [(0, 0), (1, 0)],
        'pyrene': [(0, 0), (1, 0), (0, 1), (1, 1)],
    }
    for M in range(3):
        objects['coronene' + (f'+{M}' if M else '')] = [
            (i, j) for i in range(-1-M, 2+M) for j in range(-1-M, 2+M)
            if abs(i-j) <= 1+M
        ]
    for M in [10, 20]:
        objects[f'strip-{M}'] = [
            (i, j) for i in range(M) for j in range(2)
        ]
    ctx.load_tool('aims')
    v1 = np.array((a, 0, 0))
    v2 = np.array((-a*sin(pi/6), a*cos(pi/6), 0))
    u1 = np.array((0, a/sqrt(3), 0))
    u2 = u1-v2
    bzcore = Molecule([Atom('C', m*u1 + n*u2) for m, n in [
        (0, 0), (1, 0), (2, 1), (2, 2), (1, 2), (0, 1)
    ]])
    monomers = {}
    for name, cores in objects.items():
        geom = Molecule()
        for m, n in cores:
            for atom in bzcore.shifted(m*v1 + n*v2):
                if not geom.atoms or geom.dist(atom) > 1e-10:
                    geom.atoms.append(atom)
        bm = geom.bondmatrix(1.3)
        for i in range(len(geom)):
            neigh = [j for j in np.flatnonzero(bm[i]) if i != j]
            if len(neigh) == 3:
                continue
            direc = sum(geom[neigh[j]].xyz - geom[i].xyz for j in range(2))
            direc = -direc/norm(direc)
            geom.atoms.append(Atom('H', geom[i].xyz + dCH*direc))
        monomers[name] = geom
    monomers['graphene'] = Crystal(
        [v1, v2, (0, 0, 1000)],
        [Atom('C', (0, 0, 0)), Atom('C', u1)]
    )
    dimers = {}
    monomers_cp = {}
    for name, geom in monomers.items():
        dimers[name, 'AA'] = geom + geom.shifted((0, 0, 3.35))
        dimers[name, 'AB'] = geom + geom.shifted(u1+u2+(0, 0, 3.35))
        geom_dummy = geom.copy()
        for atom in geom_dummy:
            atom.flags['dummy'] = True
            atom.prop = atom.prop.copy()
            atom.prop['mass'] *= 0.1
        monomers_cp[name, 'AA'] = geom_dummy + geom.shifted((0, 0, 3.35))
        monomers_cp[name, 'AB'] = geom_dummy + geom.shifted(u1+u2+(0, 0, 3.35))
    for name, geom in monomers.items():
        ctx(
            features='aims',
            templates='control.in geometry.in',
            xc='pbe',
            niter=100,
            geom=geom,
            relax_geometry='bfgs 0.005',
            k_grid=info.get(name, {}).get('k_grid'),
            relax_unit_cell='fixed_angles' if hasattr(geom, 'lattice') else None,
            basis='tight',
            aims='aims.master',
        ) * ctx.target(f'relax/{name}')
    methods = {
        'pbe': {'xc': 'pbe'},
        'scan': {'xc': 'dfauto scan'},
        'scan(a)': {'xc': 'dfauto scan\nxc_param scanx_c2x 0.3'},
        'scan(b)': {'xc': 'dfauto scan\nxc_param scanx_dx 2.0'},
        'scan(c)': {'xc': 'dfauto scan\nxc_param scanx_c1x 0.3\nxc_param scanx_c2x 0.3'},
        'lda': {'xc': 'pw-lda'},
        'b3lyp': {'xc': 'b3lyp'},
        'm06': {'xc': 'm06'},
        'tpss': {'xc': 'tpss'},
        'pbe+mbd': {'xc': 'pbe', 'vdw': 'many_body_dispersion'},
        'pbe+ts': {'xc': 'pbe', 'vdw': 'vdw_correction_hirshfeld'}
    }
    for basis in ['light', 'tight']:
        for method_key, method in methods.items():
            for name in monomers:
                if 'strip' in name:
                    continue
                # if basis == 'tight' and name in ['coronene+1', 'coronene+2'] \
                #         and method_key in ['b3lyp', 'm06']:
                #     continue
                vdw = {method['vdw']: ''} if 'vdw' in method else {}
                if hasattr(monomers[name], 'lattice'):
                    k_grid = info[name]['k_grid']
                    if 'many_body_dispersion' in vdw:
                        vacuum = info[name].get('vacuum')
                        if vacuum:
                            vdw['many_body_dispersion'] = {'vacuum': vacuum}
                else:
                    k_grid = None
                for system, db in [
                        ('monomers', monomers),
                        ('monomers-AA-cp', monomers_cp),
                        ('monomers-AB-cp', monomers_cp),
                        ('dimers-AA', dimers),
                        ('dimers-AB', dimers)
                ]:
                    if system.endswith('-cp') and not (basis == 'tight' and method_key == 'scan' and name != 'graphene'):
                        continue
                    if method_key == 'scan' and basis == 'tight' and name == 'graphene':
                        tier = 4
                    else:
                        tier = None
                    path = f'{system}/{name}/{method_key}/{basis}'
                    if 'dimers' in system or system.endswith('-cp'):
                        stacking = system.split('-')[1]
                        geom = db[name, stacking]
                    else:
                        geom = db[name]
                    ctx(
                        features=['aims', uncomment_tier],
                        templates='control.in geometry.in',
                        xc=method['xc'],
                        niter=100 if 'strip' in name else 1000,
                        geom=geom,
                        tier=tier,
                        **vdw,
                        xc_pre='pbe 8' if method_key == 'm06' else None,
                        k_grid=k_grid,
                        basis=basis,
                        aims='aims.master',
                    ) * ctx.target(path)


@feature('uncomment_tier')
def uncomment_tier(tsk):
    tier = tsk.consume('tier')
    if tier is None:
        return
    buffer = ''
    tier_now = None
    for l in tsk.inputs['control.in'].split('\n'):
        m = re.search(r'"(\w+) tier"', l) or re.search(r'(Further)', l)
        if m:
            tier_now = {'First': 1, 'Second': 2, 'Third': 3, 'Fourth': 4, 'Further': 5}[m.group(1)]
        m = re.search(r'#?(\s*(hydro|ionic) .*)', l)
        if m:
            l = m.group(1)
            if not (tier_now and tier_now <= tier):
                l = '#' + l
        if '####' in l:
            tier_now = None
        buffer += l + '\n'
    buffer = buffer.replace('sc_iter_limit', 'override_illconditioning .true.\nbasis_threshold 1e-4\nsc_iter_limit')
    tsk.inputs['control.in'] = buffer
