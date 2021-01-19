import pandas as pd
from pastml.annotation import preannotate_forest
from pastml.tree import read_tree

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--data', required=True, type=str)
    parser.add_argument('--tree', required=True, type=str)
    parser.add_argument('--log', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.tree, columns=['sexuality'])
    df = pd.read_csv(params.data, sep='\t', index_col=0)

    for n in tree.traverse('postorder'):
        if len(n.children) == 1:
            child = n.children[0]
            for gc in child.children:
                n.add_child(gc, dist=gc.dist + child.dist)
            n.remove_child(child)

    for n in tree.traverse('postorder'):
        if n.is_leaf():
            gender = None
            if n.name in df.index:
                gender = df.loc[n.name, 'gender']
                if not pd.isna(gender):
                    gender = {gender}
                else:
                    gender = None
        else:
            gender = {'M', 'F'}
            for c in n.children:
                gender &= getattr(c, 'gender')
        n.add_feature('gender', gender if gender else {'M', 'F'})

    transmissions = {t.up for t in tree}
    msm_het_mf, msm_het = 0, 0
    het_het, msm_msm, het_msm = 0, 0, 0
    nns = 0
    for n in tree.traverse():
        if n.is_leaf():
            continue
        msm_msm_addition = 0
        het_het_addition = 0
        if 'MSM' in getattr(n, 'sexuality'):
            msm_het += sum(1 / len(getattr(n, 'sexuality')) / len(getattr(_, 'sexuality'))
                                   for _ in n.children if 'HT' in getattr(_, 'sexuality'))
            msm_msm_addition += sum(1 / len(getattr(n, 'sexuality')) / len(getattr(_, 'sexuality'))
                                   for _ in n.children if 'MSM' in getattr(_, 'sexuality'))
            msm_het_mf += sum(1 / len(getattr(n, 'sexuality')) / len(getattr(_, 'sexuality')) / len(getattr(_, 'gender'))
                              for _ in n.children if 'HT' in getattr(_, 'sexuality') and 'F' in getattr(_, 'gender'))
        if 'HT' in getattr(n, 'sexuality'):
            het_het_addition += sum(1 / len(getattr(n, 'sexuality')) / len(getattr(_, 'sexuality'))
                                    for _ in n.children if 'HT' in getattr(_, 'sexuality'))
            het_msm += sum(1 / len(getattr(n, 'sexuality')) / len(getattr(_, 'sexuality'))
                                   for _ in n.children if 'MSM' in getattr(_, 'sexuality'))
        msm_msm += max(msm_msm_addition - 1, 0)
        het_het += max(het_het_addition - 1, 0)
        nns += len(n.children) - min(msm_msm_addition + het_het_addition, 1)

    with open(params.log, 'w+') as f:
        f.write('There are {} transmissions in the tree.\n'.format(nns))
        f.write('There are {} females in the data set.\n'
                .format(sum(1 / len(getattr(_, 'gender')) for _ in tree if 'F' in getattr(_, 'gender'))))
        f.write('There are {} MSM to HET transmissions in the data set, including {} to female HET.\n'
                .format(msm_het, msm_het_mf))
        f.write('There are {} MSM to MSM, {} HET to MSM, and {} HET to HET transmissions in the data set.\n'
                .format(msm_msm, het_msm, het_het))

