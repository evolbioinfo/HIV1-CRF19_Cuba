from collections import Counter

import pandas as pd
from pastml.annotation import preannotate_forest
from pastml.tree import read_tree

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_acr', required=True, type=str)
    parser.add_argument('--input_tree', required=True, type=str)
    parser.add_argument('--output_tab', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.input_tree)

    df = pd.read_csv(params.input_acr, index_col=0, sep='\t')
    df.index = df.index.map(str)
    df = df[[c for c in df.columns if c.startswith('RT') or c.startswith('PR')]]
    preannotate_forest(df, [tree])

    data = []

    with open(params.output_tab, 'w+') as f:
        f.write('DRM\tsensitive\tresistant\tTDR\tTDR clusters\tmax TDR cluster size\treversions\n')
        for col in sorted(df.columns):
            tips = list(tree)
            resistant_tips = [_ for _ in tips if getattr(_, col, None) == {'resistant'}]
            resistant_nodes = [_ for _ in tree.traverse() if not _.is_leaf() and getattr(_, col, None) == {'resistant'}]
            sensitive_nodes = [_ for _ in tree.traverse() if not _.is_leaf() and getattr(_, col, None) == {'sensitive'}]
            losses = [_ for _ in tree.traverse() if
                      not _.is_root() and getattr(_, col, None) == {'sensitive'} and getattr(_.up, col, None) == {
                          'resistant'}]

            def get_tdr_root(tip):
                while True:
                    if tip.is_root() or getattr(tip.up, col, None) != {'resistant'}:
                        return tip
                    tip = tip.up

            tip2root = {_: get_tdr_root(_) for _ in resistant_tips}
            tdr_counter = Counter(tip2root.values())

            n_clusters = len([_ for (_, c) in tdr_counter.items() if c > 1])
            max_tdr_cluster = tdr_counter.most_common(n=1)[0]
            n_res = len(resistant_tips)
            n_tdr = n_res - len(set(tip2root.values()))
            n_tips = len(tips)
            n_loss = len(losses)
            n_sensitive = (n_tips - n_res)
            f.write('{drm}\t{s} ({p_s:.0f}%)\t{r} ({p_r:.0f}%)\t{tdr} ({p_tdr:.0f}%)\t{clusters}\t{cluster_size}\t{loss}\n'
                    .format(drm=col.replace('_', ':'), s=n_sensitive, r=n_res, tdr=n_tdr,
                            clusters=n_clusters, loss=n_loss,
                            cluster_size=max_tdr_cluster[1] if max_tdr_cluster[1] > 1 else '',
                            p_r=100 * n_res / n_tips, p_s=100 * n_sensitive / n_tips, p_tdr=(100. * n_tdr / n_res) if n_res else 0))


