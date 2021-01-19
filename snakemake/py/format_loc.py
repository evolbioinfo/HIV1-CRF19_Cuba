import os

import pandas as pd
from pastml.file import PASTML_MARGINAL_PROBS_TAB
from pastml.models.f81_like import F81
from pastml.tree import read_tree

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--nwk', required=True, type=str)
    parser.add_argument('--loc', required=True, type=str)
    parser.add_argument('--acr', required=True, type=str)
    parser.add_argument('--ids', required=True, type=str)
    parser.add_argument('--label', required=True, type=str)
    parser.add_argument('--subtype', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.nwk)
    with open(params.ids, 'r') as f:
        ids = f.read().strip().split('\n')
    leaves = [_ for _ in tree if _.name in ids]
    mrca = leaves[0].get_common_ancestor(leaves)

    acr_df = pd.read_csv(params.acr, sep='\t', index_col=0)
    acr_df.index = acr_df.index.map(str)

    def get_value(node, col, mp_df):
        locs = acr_df.loc[node.name, col]
        print(locs)
        if isinstance(locs, str):
            locs = [locs]
        locs = [_ for _ in locs if not pd.isna(_)]
        mps = mp_df.loc[node.name, locs]
        if isinstance(mps, float):
            mps = [mps]
        return ' or '.join('{} ({:.2f})'.format(l, p) for (l, p) in zip(locs, mps))

    values = []
    columns = ['country', 'intregion']
    for col in columns:
        mp_df = pd.read_csv(os.path.join(os.path.dirname(params.acr),
                                         PASTML_MARGINAL_PROBS_TAB.format(model=F81, state=col)), sep='\t', index_col=0)
        mp_df.index = mp_df.index.map(str)
        values.append(get_value(mrca, col, mp_df))
        values.append(get_value(mrca.up, col, mp_df))

    with open(params.loc, 'w') as f:
        f.write('data set\tsubtype\t{}\n{}\t{}\t{}'
                .format('\t'.join(['{c}: CRF_19 MRCA\t{c}: non-CRF_19 MRCA'.format(c=c) for c in columns]),
                        params.label, params.subtype, '\t'.join(values)))
