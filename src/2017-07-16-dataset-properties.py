import sys
sys.path.append('../calc/2017-01-23-all-vdw-sets-3')
import caf  # noqa
import vdwsets.vdwsets as vdw
import pandas as pd

# ::>


def get_dist(m1, m2):
    return m1.dist(m2)


(
    pd.DataFrame([
        (dsname, *k, get_dist(*ds.geoms[c.fragments['complex']].get_fragments(1.1)))
        for dsname, ds in vdw.get_all_datasets(include='S66x8 S12L'.split()).items()
        for k, c in ds.clusters.items()
    ], columns='dataset system dist real_dist'.split())
    .loc[lambda df: [df.real_dist.idxmax(), df.real_dist.idxmin()]]
)

# ::>


def get_ref(refs):
    return refs.get('ref') or refs.get('ref-alberto') or refs.get('ref-grimme')


(
    pd.DataFrame([
        (dsname, system, dist[0] if dist else 1., get_ref(c.energies))
        for dsname, ds in vdw.get_all_datasets(include='S22 S66x8 S12L'.split()).items()
        for (system, *dist), c in ds.clusters.items()
    ], columns='dataset system dist ene'.split())
    .query('dist == 1')
    .loc[lambda df: [df.ene.idxmax(), df.ene.idxmin()]]
)
