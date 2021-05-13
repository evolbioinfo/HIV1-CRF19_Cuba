import pandas as pd
from ete3 import Tree
from ete3.parser.newick import write_newick


def _get_date(d, exact_only=False):
    if pd.notnull(d):
        first_jan_this_year = pd.datetime(year=d.year, month=1, day=1)
        day_of_this_year = d - first_jan_this_year
        first_jan_next_year = pd.datetime(year=d.year + 1, month=1, day=1)
        days_in_this_year = first_jan_next_year - first_jan_this_year
        date = d.year + day_of_this_year / days_in_this_year
        if exact_only and (date == d.year or d.day == 1):
            return None
        return date
    else:
        return None


def read_tree(nwk):
    tree = None
    for format in range(10):
        try:
            tree = Tree(nwk, format=format)
            break
        except:
            continue
    return tree


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--data', required=True, type=str)
    parser.add_argument('--in_tree', required=True, type=str)
    parser.add_argument('--out_tree', required=True, type=str)
    parser.add_argument('--date_col_tip', required=True, type=str)
    parser.add_argument('--date_col_cut', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.in_tree)

    date_df = pd.read_table(params.data, index_col=0)[[params.date_col_tip, params.date_col_cut]]
    date_df.index = date_df.index.map(str)
    date_df = date_df.loc[[_.name for _ in tree], :]

    date_df[params.date_col_tip] = pd.to_datetime(date_df[params.date_col_tip].astype(str).str.replace('.0', '', regex=False),
                                         infer_datetime_format=True, errors='coerce').map(_get_date)
    date_df[params.date_col_cut] = pd.to_datetime(date_df[params.date_col_cut].astype(str).str.replace('.0', '', regex=False),
                                         infer_datetime_format=True, errors='coerce').map(_get_date)

    tip = next(_ for _ in tree if not pd.isna(date_df.loc[_.name, params.date_col_tip]))
    tip.add_feature('date', date_df.loc[tip.name, params.date_col_tip])
    while not tip.is_root():
        tip.up.add_feature('date', getattr(tip, 'date') - tip.dist)
        tip = tip.up
    print('Root date of tree {} is {}'.format(params.in_tree, getattr(tree, 'date')))
    for _ in tree.traverse('preorder'):
        if not _.is_root():
            _.add_feature('date', getattr(_.up, 'date') + _.dist)

    for _ in tree:
        cut_date = date_df.loc[_.name, params.date_col_cut]
        tip_date = getattr(_, 'date')
        if not pd.isna(cut_date) and cut_date < tip_date:
            delta = tip_date - cut_date
            if 0 <= delta <= _.dist:
                _.dist -= delta
            else:
                node = _
                while True:
                    delta -= node.dist
                    parent = node.up
                    parent.remove_child(node)
                    grandparent = parent.up
                    for sibling in list(parent.children):
                        sibling.dist += parent.dist
                        grandparent.add_child(sibling)

                    if parent.dist >= delta:
                        grandparent.remove_child(parent)
                        grandparent.add_child(_, dist=parent.dist - delta)
                        break
                    else:
                        node = parent
    nwk = write_newick(tree, format_root_node=True)
    with open(params.out_tree, 'w+') as f:
        f.write('%s\n' % nwk)


