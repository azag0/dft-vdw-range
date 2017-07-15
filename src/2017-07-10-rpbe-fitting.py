# ::hide
import os
from itertools import product
from collections import OrderedDict
import numpy as np
import pandas as pd
from math import ceil
from scipy.interpolate import interp1d
from scipy.optimize import minimize, brentq
# ::>


# ::hide
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
plt.style.use('seaborn-whitegrid')
mpl.rc('font', family='serif', serif='STIXGeneral')
mpl.rc('mathtext', fontset='stix')
import seaborn as sns
sns.set(context='paper', style='whitegrid', font='serif', rc={
    'axes.formatter.useoffset': False,
    'font.serif': 'STIXGeneral',
    'mathtext.fontset': 'stix'
})
# ::%matplotlib inline
# ::%config InlineBackend.figure_format = 'svg'
# ::>


# ::hide
data = {}
with pd.HDFStore('../data/all-data.h5') as store:
    for key in 'scf mbd dftd3_bj vv10 flakes bz_tests'.split():
        data[key] = store[key]

# filter out data points where distance = 2 and reference energy is nonnegative
mask = (data['scf'].index.get_level_values('dist') == 2.) & (data['scf']['ref'] >= 0)
data['scf'].drop(data['scf'].loc[mask].index, inplace=True)
# ::>


# ::hide
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


data['scf'] = data['scf'] \
    .assign(base=lambda x: abs(x['ref'])) \
    .groupby(level='name system xc'.split()) \
    .apply(mod_ref1)
# ::>


# ::hide
def savefig(name, **kwargs):
    sns.plt.savefig(
        os.path.expanduser(f'../media/{name}.pdf'),
        transparent=True,
        **{'bbox_inches': 'tight', **kwargs}
    )


def get_stat(df):
    return pd.Series({
        'STD': df['reldelta'].std(),
        'mean': df['reldelta'].mean(),
        'MARE': abs(df['reldelta']).mean(),
        'median': df['reldelta'].median(),
        'MAE': abs(df['delta']).mean(),
    })


def normalize_names(col):
    names = {
        'scan': 'SCAN',
        'scan0': 'SCAN0',
        'scan(a)': "SCAN*\n($c_{2\\mathrm{x}}=0.3$)",
        'scan(b)': "SCAN*\n($d_\\mathrm{x}=2$)",
        'b3lyp': 'B3LYP',
        'lda': 'LDA',
        'm06': 'M06',
        'pbe': 'PBE',
        'pbe0': 'PBE0',
        'tpss': 'TPSS',
        'pbe+mbd': 'PBE+MBD',
    }
    if col.name.lower() == 'xc':
        return [names.get(k, k) for k in col]
    return col


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


# ::>
# **Relative error distributions on S66x8.** Obviously, RPBE is not only much
# more shortranged than PBE, but even than B3LYP. You cannot go really much
# lower with $\beta$ than 0.65, so the RPBE case shows the best possible
# scenario, and you see that it's extremely repulsive around equilibrium.


# ::hide
def get_name(row):
    name = row['xc'].upper()
    name += '\n({})'.format(f'$\\beta={row["beta"]}$' if row['beta'] != 10. else 'no vdW')
    return name


keys = [
    ('lda', 10.), ('pbe', .83), ('pbe0', .83),
    ('rpbe', .67),
    ('b3lyp', .67), ('scan', 1.1), ('scan', 1.4),
]
with sns.color_palette(list(reversed(sns.color_palette('coolwarm', 8)))):
    g = sns.factorplot(
        data=pd.concat([merge_scf_vdw(
            'mbd', ('S66x8', slice(None), slice(None)), (xc,), (6, beta), 'a beta'.split()
        ).reset_index() for xc, beta in keys]).assign(
            xc_beta=lambda x: x.apply(get_name, axis=1)
        ),
        kind='box',
        x='xc_beta',
        y='reldelta',
        hue='dist',
        whis=2.5,
        aspect=3.2,
        size=2.,
        margin_titles=True
    )
sns.plt.ylim(-.5, .55)
g.set_xlabels('')
g.set_ylabels(r'$\Delta E_i/E_i^\mathrm{ref}$')
g.set(yticks=[-.3, -.1, 0, .1, .3])
g.set_yticklabels(['$-30\%$', '$-10\%$', '0%', '10%', '30%'])
# ::>


# ::hide
def interp_vv10(bs):
    def func(r):
        f = interp1d(r.index, r)
        return pd.Series({b: ene for b, ene in zip(bs, f(bs))})
    return data['vv10']['ene'] \
        .unstack().apply(lambda r: r['nlc'] if r.name[0] == 'S12L' else r['vdw']-r['base'], axis=1) \
        .unstack().apply(func, axis=1).rename_axis('b', axis=1).stack().to_frame('ene')


# ::>
# **Fitting damping parameters on S66.** This just confirms the picture above:
# RPBE is extremely shortranged. Neither MBD nor D3 can handle that.
#
# But interestingly, VV10 can handle that with extremely small damping
# parameter. Now the question is how transfarable that is because at this point
# the nonlocal functional is covering midrange parts of the interaction for
# which it really isn't designed.
#
# I don't have TS evaluted, but I'm pretty sure it behaves very similarly to MBD
# in this regard, so I wouldn't expect RPBE+TS would be any good.


# ::hide
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
            ax.set_ylim((-0.4, 0.4))
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
    funcs = ['pbe', 'scan', 'rpbe']
    setnames = ['S66x8']
    methods = OrderedDict((k, v) for k, *v in [
        ('mbd', 'beta', r'$\beta^\mathrm{MBD}$', [.8, 1, 1.2]),
        ('vv10', 'b', r'$b^\mathrm{VV10}$', [10, 20]),
        ('d3', 'a2', r'$a_2^\mathrm{D3}$', [2, 8]),
    ])
    fig = plt.figure(figsize=(4, 4))
    axes = get_axes(
        fig, funcs=funcs, setnames=setnames, methods=methods, ncols=1,
        hspace=0.14, wspace=0.08, top_all=0.9
    )
    my_data = get_data(funcs=funcs, setnames=setnames)
    handles = {}
    for method, df in my_data.items():
        for (func, setname), gr in df.groupby(['xc', 'name']):
            plot_stat2(axes[func, method, setname], gr, methods[method][0], handles)
    fig.legend(handles.values(), handles.keys(), 'upper center', ncol=3, frameon=True)
    return fig


plot_all_funcs()
# ::>
