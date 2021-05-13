import pandas as pd
from pastml.acr import annotate_dates
from pastml.annotation import preannotate_forest
from pastml.tree import read_tree, collapse_zero_branches

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--acr', required=True, type=str)
    parser.add_argument('--input_tree', required=True, type=str)
    parser.add_argument('--output_tree', required=True, type=str)
    parser.add_argument('--output_log', required=True, type=str)
    parser.add_argument('--value', required=True, type=str)
    parser.add_argument('--root_date', required=True, type=float)
    parser.add_argument('--ids', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.input_tree)
    collapse_zero_branches([tree])
    annotate_dates([tree], root_dates=[params.root_date])
    tab_df = pd.read_csv(params.acr, sep='\t', header=0, index_col=0)
    col = preannotate_forest(tab_df, [tree])[0]
    for n in tree.traverse('levelorder'):
        if not n.is_leaf() and getattr(n, col) == {params.value}:
            n.write(outfile=params.output_tree, format=3)
            with open(params.output_log, 'w') as f:
                f.write('{}'.format(getattr(n, 'date')))
            with open(params.ids, 'w') as f:
                f.write('\n'.join(_.name for _ in n))
            break

