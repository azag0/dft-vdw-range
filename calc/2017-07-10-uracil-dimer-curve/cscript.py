#!/usr/bin/env python3
try:
    import pymbd
except:
    pass
import numpy as np
import pandas as pd

import caf
from caflib.Configure import function_task
from caflib.Tools import geomlib2
from caflib.Tools.aims import AimsTask

ev = 27.2107
kcal = 627.503

paths = ['<>/<>/*/<><,/vols>']

base = dict(
    spin='none',
    relativistic='none',
    occupation_type=('gaussian', 0.02),
    mixer='pulay',
    n_max_pulay=8,
    charge_mix_param=0.3,
    sc_accuracy_eev=1e-3,
    sc_accuracy_rho=1e-5,
    sc_accuracy_etot=1e-6,
    sc_iter_limit=100,
)
xcs = [
    ('lda', 'pw-lda'),
    ('blyp', 'blyp'),
    ('pbe', 'pbe'),
    ('rpbe', 'rpbe'),
    ('scan', 'dfauto scan'),
    ('tpss', 'tpss'),
    ('m06-l', 'm06-l'),
]
systems = [
    ('uracil-dimer', 'tight', [-.6, -.5, -.4, -.2, 0, .2, .4, .7, 1, 1.4, 1.8, 2.1, 2.5, 3., 3.5, 5, 100]),
    ('bucky-catcher', 'tight', [-.4, -.2, 0, .2, .4, .7, 1, 1.4, 1.8, 2.5, 3.5, 5, 100]),
]


def run(ctx):
    data = []
    volumes = {}
    eq_dists = {}
    for system, basis, dists in systems:
        mon1, mon2 = geomlib2.readfile(f'{system}.xyz').get_fragments()
        direc = mon2.cms-mon1.cms
        eq_dist = eq_dists.setdefault(system, np.linalg.norm(direc))
        direc = direc/eq_dist
        for xc_label, xc in xcs:
            for d in dists:
                geoms = [
                    ('complex', mon1 + mon2.shifted(d*direc)),
                    ('fragment-1', mon1 + mon2.shifted(d*direc).ghost()),
                    ('fragment-2', mon1.ghost() + mon2.shifted(d*direc)),
                ]
                for geom_label, geom in geoms:
                    dft = ctx(
                        klass=AimsTask,
                        geom=geom,
                        aims='aims.b4baa15',
                        basis=basis,
                        tags=dict(
                            **base,
                            xc=xc,
                            output='hirshfeld' if xc_label == 'pbe' else None,
                        )
                    )
                    calc = get_energy(
                        ('calc', dft.outputs['run.out']),
                        target=f'{system}/{xc_label}/{d}/{geom_label}',
                        ctx=ctx
                    )
                    if calc.finished:
                        data.append((system, xc_label, eq_dist+d, geom_label, calc.result/ev*kcal))
                    if xc == 'pbe':
                        calc = get_volumes(
                            ('calc', dft.outputs['run.out']),
                            target=f'{system}/{xc_label}/{d}/{geom_label}/vols',
                            ctx=ctx
                        )
                        if calc.finished:
                            volumes[system, geom_label, eq_dist+d] = geom, calc.result
    return (
        pd.DataFrame(data, columns='system xc d geom ene'.split()),
        volumes,
        eq_dists
    )


@function_task
def get_energy(output):
    with open(output) as f:
        next(l for l in f if l == '  Self-consistency cycle converged.\n')
        return float(next(l for l in f if 'Total energy uncorrected' in l).split()[5])


@function_task
def get_volumes(output):
    from math import nan

    ratios = []
    with open(output) as f:
        next(l for l in f if 'Performing Hirshfeld analysis' in l)
        for line in f:
            if not line.strip():
                break
            if 'Free atom volume' in line:
                free = float(line.split()[-1]) or nan
            elif 'Hirshfeld volume' in line:
                hirsh = float(line.split()[-1])
                ratios.append(hirsh/free)
    return ratios


# ::>
# <!-- END CSCRIPT -->

import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns  # noqa
matplotlib.pyplot.style.use('seaborn-whitegrid')
matplotlib.rc('font', family='serif', serif='STIXGeneral')
matplotlib.rc('mathtext', fontset='stix')
# ::%matplotlib inline
# ::%config InlineBackend.figure_format = 'svg'
# ::%config InlineBackend.print_figure_kwargs = {'bbox_inches': 'tight', 'dpi': 300}

pymbd.lib.init_grid(15)


def savefig(fig, name):
    fig.savefig(
        f'../../media/{name}.pdf',
        bbox_inches='tight',
        dpi=600,
    )


# ::>

results, volumes, eq_dists = run(caf.get_context())

# ::>

mbd_data = []
for (system, fragment, d), (geom, vols) in volumes.items():
    vols = np.array(vols)
    vols = vols[~np.isnan(vols)]
    ene = pymbd.mbd_rsscs(geom.xyz/pymbd.bohr, geom.species, vols, 0.83)
    mbd_data.append((system, fragment, 'pbe', d, ene*kcal))
mbd_data = pd.DataFrame(mbd_data, columns='system fragment xc d ene'.split())

# ::>


def get_ene_int(df):
    return df.ene['complex']-df.ene['fragment-1']-df.ene['fragment-2']


ene_int = results.set_index('system xc d geom'.split()).unstack().apply(get_ene_int, axis=1)
ene_int = pd.concat([
    ene_int,
    ene_int.to_frame('ene').reset_index().join(
        mbd_data.set_index('system fragment xc d'.split())
        .unstack('fragment').apply(get_ene_int, axis=1).to_frame('ene'),
        on='system xc d'.split(),
        rsuffix='_mbd',
        how='inner'
    )
    .assign(
        ene=lambda df: df.ene+df.ene_mbd,
        xc=lambda df: df.xc + '+mbd',
    ).drop('ene_mbd', axis=1).set_index('system xc d'.split()).ene
])

# ::>

from mpl_toolkits.axes_grid.inset_locator import (
    zoomed_inset_axes, mark_inset, InsetPosition
)
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from itertools import product


dfs = ene_int.to_frame('ene').dropna().reset_index().query('d < 100').groupby('system xc'.split())
cmap = iter(sns.husl_palette(9))
xc_labels = 'BLYP RPBE TPSS PBE M06-L SCAN LDA PBE+MBD'.split()
colors = {xc: next(cmap) for xc in 'B3LYP TPSS PBE SCAN LDA M06'.split()}
fig, axes = plt.subplots(2, 1, figsize=(3.7, 4), gridspec_kw=dict(top=.8))
axes_in = [zoomed_inset_axes(axes[i], 1) for i in range(2)]
handles = {}
for ((system, *_), (ax0, ax1)), xc in product(zip(systems, zip(axes, axes_in)), xc_labels):
    df = dfs.get_group((system, xc.lower()))
    minmax = df.d.min()+1e-3, df.d.max()-1e-3
    fene = interp1d(df.d, df.ene, kind='cubic')
    d0 = minimize(fene, 7, bounds=[minmax]).x
    d = np.linspace(*minmax, 200)
    ene = fene(d)
    col_label = xc.split('+', 1)[0]
    col = colors.get(col_label) or colors.setdefault(col_label, next(cmap))
    handles[xc], = ax0.plot(d, ene, color=col, lw=1, ls='--' if '+' in xc else None)
    # ax0.plot(df.d, df.ene, ls='none', marker='o', color=col, fillstyle='none', mew=1, ms=3)
    if ax1:
        ax1.plot(d, ene, color=col, lw=1, ls='--' if '+' in xc else None)
    if d0 < 7:
        ax0.plot(d0, fene(d0), marker='o', color=col, fillstyle='none', mew=1.3, ms=5)
    try:
        ax0.plot(eq_dists[system], dfs.get_group((system, xc.lower() + '+')).ene.iloc[0], marker='s', color=col, fillstyle='none', mew=1.3, ms=5)
    except:
        pass

for i in range(2):
    axes[i].set_ylabel('$E$ [kcal/mol]')
axes[1].set_xlabel('center-of-mass distance [$\mathrm{\\AA}$]')
axes[0].set_xlim(2.8, 5.7)
axes[0].set_ylim(-12, 0)
axes[0].set_xticks([eq_dists['uracil-dimer'], 5])
axes[0].set_xticklabels([np.round(eq_dists['uracil-dimer'], 1), 5])
axes[0].set_yticks([0, -9.8])
axes[0].set_yticklabels([0, '$-9.8$'])
axes[0].plot(eq_dists['uracil-dimer'], -9.8, marker='o', color='black', fillstyle='none', mew=2, ms=5)
axes_in[0].set_xlim(4.9, 5.4)
axes_in[0].set_ylim(-2.1, -.5)
axes_in[0].set_xticks([])
axes_in[0].set_yticks([])
axes_in[0].set_axes_locator(InsetPosition(axes[0], (.6, .04, .37, .6)))
mark_inset(axes[0], axes_in[0], 2, 1)
imbox = OffsetImage(plt.imread('../../media/uracil-dimer.png'), zoom=.065)
axes[0].add_artist(AnnotationBbox(imbox, xy=(4.15, -8.7), frameon=False))
axes[1].set_xlim(4.5, 10)
axes[1].set_ylim(-35, 5)
axes[1].set_xticks([eq_dists['bucky-catcher'], 9])
axes[1].set_xticklabels([np.round(eq_dists['bucky-catcher'], 1), 9])
axes[1].set_yticks([-30, 0])
axes[1].plot(eq_dists['bucky-catcher'], -30, marker='o', color='black', fillstyle='none', mew=2, ms=5)
axes_in[1].set_xlim(8.1, 8.9)
axes_in[1].set_ylim(-3, 2)
axes_in[1].set_xticks([])
axes_in[1].set_yticks([0])
axes_in[1].set_axes_locator(InsetPosition(axes[1], (.6, .04, .37, .6)))
mark_inset(axes[1], axes_in[1], 2, 1)
imbox = OffsetImage(plt.imread('../../media/bucky-catcher.png'), zoom=.065)
axes[1].add_artist(AnnotationBbox(imbox, xy=(6.8, -25), frameon=False))
fig.legend(handles.values(), handles.keys(), 'upper center', ncol=3, frameon=True)
savefig(fig, 'range-curves')
