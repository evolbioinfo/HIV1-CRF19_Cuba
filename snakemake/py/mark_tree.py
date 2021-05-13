import pandas as pd
from ete3 import Tree, TreeNode
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

    parser.add_argument('--data', default='/home/azhukova/atlas/projects/cuba_atlas/data/datasets/D_CRF_19/metadata.tab', type=str)
    parser.add_argument('--in_tree', default='/home/azhukova/atlas/projects/cuba_atlas/data/datasets/D_CRF_19/timetree.lsd2.collapsed.nwk', type=str)
    parser.add_argument('--out_tree', default='/home/azhukova/atlas/projects/cuba_atlas/data/datasets/D_CRF_19/timetree.lsd2.collapsed.diag.nwk', type=str)
    parser.add_argument('--date_col_tip', default='sample_date', type=str)
    parser.add_argument('--date_col_cut', default='diagnostics_date', type=str)
    params = parser.parse_args()

    tree = read_tree(params.in_tree)
    date_df = pd.read_csv(params.data, index_col=0, sep='\t')

    if params.date_col_cut in date_df.columns:
        date_df.index = date_df.index.map(str)
        date_df = date_df.loc[[_.name for _ in tree], :]

        date_df[params.date_col_tip] = pd.to_datetime(
            date_df[params.date_col_tip].astype(str).str.replace('.0', '', regex=False),
            infer_datetime_format=True, errors='coerce').map(_get_date)
        date_df[params.date_col_cut] = pd.to_datetime(
            date_df[params.date_col_cut].astype(str).str.replace('.0', '', regex=False),
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

        for tip in tree:
            cut_date = date_df.loc[tip.name, params.date_col_cut]
            tip_date = getattr(tip, 'date')
            if not pd.isna(cut_date) and cut_date < tip_date:
                node = tip
                while getattr(node, 'date') - node.dist > cut_date:
                    node = node.up
                delta = getattr(node, 'date') - cut_date
                mark_node = TreeNode(name=tip.name + '_diagnostics')
                mark_node.add_feature('date', cut_date)
                parent = node.up
                parent.remove_child(node)
                parent.add_child(mark_node, dist=node.dist - delta)
                mark_node.add_child(node, dist=delta)

    # if there is a branch with several 'diagnostics' nodes on it -- fail,
    # as it means the dating went wrong, and the divergence date is too recent
    todo = [tree]
    while todo:
        n = todo.pop()
        fake_nodes = []
        while len(n.children) == 1:
            fake_nodes.append(n)
            n = n.children[0]
        if len(fake_nodes) > 1:
            raise ValueError('The diagnostics dates for {} are in conflict with the tree structire '
                             'that suggest that at that time only one of those individuals was infected.'
                             .format(', '.join(fake_node.name.replace('_diagnostics', '') for fake_node in fake_nodes)))
            # # potential fix that simply prunes the conflicting tips
            # fake_node = fake_nodes[0]
            # tip = next(_ for _ in fake_node if _.name == fake_node.name.replace('_diagnostics', ''))
            # to_remove = [_ for _ in fake_node if _ != tip]
            # for node in to_remove:
            #     parent = node.up
            #     parent.remove_child(node)
            #     # If the parent node has only one child now, merge them.
            #     if len(parent.children) == 1 and parent != fake_node:
            #         brother = parent.children[0]
            #         brother.dist += parent.dist
            #         grandparent = parent.up
            #         grandparent.remove_child(parent)
            #         grandparent.add_child(brother)
            #     else:
            #         while parent.is_leaf():
            #             grandparent = parent.up
            #             grandparent.remove_child(parent)
            #             parent = grandparent
            # while tip.up != fake_node:
            #     parent = tip.up
            #     tip.dist += parent.dist
            #     grandparent = parent.up
            #     grandparent.remove_child(parent)
            #     grandparent.add_child(tip)

        todo.extend(n.children)

    nwk = write_newick(tree, format_root_node=True)
    with open(params.out_tree, 'w+') as f:
        f.write('%s\n' % nwk)
