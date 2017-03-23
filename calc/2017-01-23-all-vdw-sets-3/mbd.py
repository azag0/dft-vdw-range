import json
import numpy as np
from itertools import product


def run(ctx, mbd):
    mbd.init_grid(15)
    mode = ''
    with open('input.json') as f:
        inp = json.load(f)
    coords = np.array(inp['coords'])/mbd.bohr
    species = inp['species']
    if 'lattice' in inp:
        lattice = np.array(inp['lattice'])/mbd.bohr
        mode += 'C'
        k_grid = mbd.make_k_grid(mbd.make_g_grid(*inp['k_grid']), lattice)
    else:
        lattice = None
    with open('dft.json') as f:
        volumes = json.load(f)['hirsh']
    alpha_0, C6, R_vdw = ctx.scale_hirsh(volumes, *ctx.get_free_atom_data(species))
    results = []
    betas = [
        .57, .6, .63, .67, .7, .73, .77, .8, .83, .87, .9, .95, 1.,
        1.1, 1.2, 1.23, 1.27, 1.3, 1.4, 1.6, 1.8, 2., 3.
    ]
    a_s = [5, 6, 7, 9, 11, 13, 15]
    for beta, a in product(betas, a_s):
        alpha_scs_dyn = mbd.run_scs(
            mode, 'fermi,dip,gg', coords,
            mbd.alpha_dynamic_ts_all('C', mbd.n_grid_omega, alpha_0, c6=C6),
            unit_cell=lattice,
            r_vdw=R_vdw, beta=beta, a=a
        )
        C6_scs = mbd.get_c6_from_alpha(alpha_scs_dyn)
        R_vdw_scs = R_vdw*(alpha_scs_dyn[0]/alpha_0)**(1/3)
        if 'lattice' in inp:
            ene = mbd.get_reciprocal_mbd_energy(
                mode.replace('C', 'R'),
                'fermi,dip',
                coords,
                alpha_scs_dyn[0],
                mbd.omega_eff(C6_scs, alpha_scs_dyn[0]),
                k_grid,
                unit_cell=lattice,
                r_vdw=R_vdw_scs,
                beta=beta,
                a=a
            )[0]
        else:
            ene = mbd.get_single_mbd_energy(
                mode,
                'fermi,dip',
                coords,
                alpha_scs_dyn[0],
                mbd.omega_eff(C6_scs, alpha_scs_dyn[0]),
                r_vdw=R_vdw_scs,
                beta=beta,
                a=a
            )[0]
        results.append((beta, a, ene))
    if ctx.myid == 0:
        with open('mbd.json', 'w') as f:
            json.dump(results, f)
