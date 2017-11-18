from itertools import product
from collections import OrderedDict
from math import ceil

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize, brentq

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
plt.style.use('seaborn-whitegrid')
mpl.rc('font', family='serif', serif='STIXGeneral')
mpl.rc('mathtext', fontset='stix')
# ::%matplotlib inline
# ::%config InlineBackend.figure_format = 'svg'
# ::>


data = {}
with pd.HDFStore('../data/all-data.h5') as store:
    for key in 'scf mbd dftd3_bj vv10 flakes bz_tests'.split():
        data[key] = store[key]
# ::>


# set base energy (used to get relative energy) to abs(ref)
# and flatten the dip around ene_ref = 0
def mod_ref1(df):
    if len(df) == 1:
        return df
    ibase = df['ref'].idxmin()
    base = df.loc[ibase, 'base']
    dist = ibase[2]
    mask = (df.index.get_level_values('dist') < dist) & (df['base'] < base)
    df.loc[mask, 'base'] = base
    return df


data['scf'] = (
    data['scf']
    .assign(base=lambda x: abs(x['ref']))
    .groupby('name system xc'.split())
    .apply(mod_ref1)
)
# ::>


def savefig(fig, name):
    fig.savefig(
        f'../media/{name}.pdf',
        transparent=True,
        bbox_inches='tight',
    )


def get_stat(df):
    return pd.Series({
        'STD': df['reldelta'].std(),
        'mean': df['reldelta'].mean(),
        'MARE': abs(df['reldelta']).mean(),
        'median': df['reldelta'].median(),
        'MAE': abs(df['delta']).mean(),
    })


def merge_scfs_vdws(vdw, baseidx, scfidxs, vdwidxs, params):
    data_vdw = data[vdw] if type(vdw) is str else vdw
    df = pd.merge(
        pd.concat(data['scf'].loc[baseidx + idx, :] for idx in scfidxs).reset_index(),
        pd.concat(data_vdw.loc[baseidx + idx, :] for idx in vdwidxs).reset_index(),
        on='name system dist'.split(),
        how='inner',
        suffixes=('_scf', '_vdw')
    )
    df.set_index('name xc system dist'.split() + params, inplace=True)
    df['ene'] = df['ene_scf']+df['ene_vdw']
    df['delta'] = df['ene']-df['ref']
    df['reldelta'] = df['delta']/df['base']
    return df


def merge_scf_vdw(vdw, baseidx, scfidx, vdwidx, params):
    return merge_scfs_vdws(vdw, baseidx, [scfidx], [vdwidx], params)


def interp_vv10(bs):
    def func(r):
        f = interp1d(r.index, r)
        return pd.Series({b: ene for b, ene in zip(bs, f(bs))})
    return data['vv10']['ene'] \
        .unstack().apply(lambda r: r['nlc'] if r.name[0] == 'S12L' else r['vdw']-r['base'], axis=1) \
        .unstack().apply(func, axis=1).rename_axis('b', axis=1).stack().to_frame('ene')
# ::>


# ::%%time
def get_data(*, funcs, setnames):
    idxs = (setnames, slice(None), 1.), (funcs,)
    return {
        'mbd': merge_scf_vdw('mbd', *idxs, (6, slice(.7, 1.5)), ['a', 'beta']),
        'd3': merge_scf_vdw('dftd3_bj', *idxs, (.55, slice(2., 11.), 0.), ['a1', 'a2', 's8']),
        'vv10': merge_scf_vdw(interp_vv10(np.linspace(4.5, 22., 10)), *idxs, (), ['b'])
    }


def get_axes(fig, *, funcs, setnames, methods, ncols, hspace, wspace, top_all):
    axes = {}
    nsetnames = len(setnames)
    nmethods = len(methods)
    nrows = ceil(len(funcs)/ncols)
    height = (top_all-(nrows-1)*hspace)/nrows
    width = (1-(ncols-1)*wspace)/ncols
    drow = height+hspace
    dcol = width+wspace
    for func, (irow, icol) in zip(funcs, product(range(nrows), range(ncols))):
        top, left = top_all-irow*drow, 0+icol*dcol
        bottom, right = top-height, left+width
        gs = GridSpec(
            nsetnames, nmethods, top=top, bottom=bottom, left=left,
            right=right, hspace=0.05, wspace=0.05)
        fig.text(left+width/2, top+0.01, func.upper(), ha='center')
        for (i, setname), (j, (method, (param, param_label, xticks))) in \
                product(enumerate(setnames), enumerate(methods.items())):
            ax = fig.add_subplot(gs[i, j])
            ax.set_ylim((-0.3, 0.3))
            ax.set_xticks(xticks)
            if i == nsetnames-1:
                ax.set_xlabel(param_label)
            else:
                for tl in ax.xaxis.get_ticklabels():
                    tl.set_visible(False)
            ax.set_yticks([-.1, 0, .1])
            if j == 0:
                ax.set_yticklabels(['$-10\%$', '0%', '10%'])
            else:
                for tl in ax.yaxis.get_ticklabels():
                    tl.set_visible(False)
            axes[func, method, setname] = ax
    return axes


def plot_stat2(ax, df, param, handles):
    df = df.groupby(param).apply(get_stat).reset_index()
    mean = interp1d(df[param], df['mean'], kind='cubic')
    std = interp1d(df[param], df['STD'], kind='cubic')
    param_vals = np.linspace(df[param].min(), df[param].max())
    line, = ax.plot(param_vals, mean(param_vals), '--', color='black')
    handles['MRE'] = line
    line, = ax.plot(param_vals, std(param_vals), '-', color='black')
    handles['SDRE'] = line
    bounds = [(df[param].min()+1e-6, df[param].max()-1e-6)]
    opt_std = minimize(std, np.mean(bounds), bounds=bounds).x[0]
    if np.sign(mean(bounds[0][0])) != np.sign(mean(bounds[0][1])):
        opt_mean = brentq(mean, *bounds[0])
        ax.axvline(opt_mean, color='black', ls='dotted', lw=1)
    ax.axhline(color='gray', lw=1.5)
    line = ax.axvline(opt_std, color='black', ls='dotted', lw=1)
    handles['optimal MRE/SDRE'] = line


def plot_all_funcs():
    funcs = ['lda', 'pbe', 'pbe0', 'b3lyp', 'scan', 'm06-l']
    setnames = ['S66x8', 'X23', 'S12L']
    methods = OrderedDict((k, v) for k, *v in [
        ('mbd', 'beta', r'$\beta^\mathrm{MBD}$', [.8, 1, 1.2]),
        ('vv10', 'b', r'$b^\mathrm{VV10}$', [10, 20]),
        ('d3', 'a2', r'$a_2^\mathrm{D3}$', [2, 8]),
    ])
    fig = plt.figure(figsize=(7, 9))
    axes = get_axes(
        fig, funcs=funcs, setnames=setnames, methods=methods, ncols=2,
        hspace=0.08, wspace=0.08, top_all=0.92
    )
    my_data = get_data(funcs=funcs, setnames=setnames)
    handles = {}
    for method, df in my_data.items():
        for (func, setname), gr in df.groupby(['xc', 'name']):
            plot_stat2(axes[func, method, setname], gr, methods[method][0], handles)
    fig.legend(handles.values(), handles.keys(), 'upper center', ncol=3, frameon=True)
    return fig


savefig(plot_all_funcs(), 'param-fitting-all')
# ::>
