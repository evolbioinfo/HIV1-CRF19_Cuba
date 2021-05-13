import re
import os

from Bio import Phylo

DATE_REGEX = r'[0-9]{4}[.0-9]*'
RATE_REGEX = r'rate\s*([\.\de+-]+\s\[[\.\de+-]+;\s*[\.\de+-]+\])'

DATE_CI_FORMAT = '{:.2f} [{:.2f} - {:.2f}]'


def format_date_ci(confidence):
    return DATE_CI_FORMAT.format(*[float(_) for _ in re.findall(DATE_REGEX, confidence)])


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--nexus', required=True, type=str)
    parser.add_argument('--log', required=True, type=str)
    parser.add_argument('--dates', required=True, type=str)
    parser.add_argument('--ids', required=False, type=str, default=None)
    parser.add_argument('--label', required=True, type=str)
    parser.add_argument('--subtype', required=True, type=str)
    params = parser.parse_args()

    with open(params.log, 'r') as f:
        for line in f.readlines():
            rate = re.findall(RATE_REGEX, line)
            if rate:
                rate = rate[0]
                break

    with open(params.nexus, 'r') as f:
        tree = f.read()
    tree = re.sub(r'CI_date="({})\(({}),({})\)"'.format(DATE_REGEX, DATE_REGEX, DATE_REGEX), r'CI_date="\1 \2 \3"', tree)
    temp = params.nexus + '.temp'
    with open(temp, 'w') as f:
        f.write(tree)
    tree = Phylo.read(temp, 'nexus')
    os.remove(temp)
    if params.ids is None:
        mrca = tree.root
        mrca_ci = format_date_ci(tree.root.branch_length)
    else:
        with open(params.ids, 'r') as f:
            ids = f.read().strip().split('\n')
        mrca = tree.common_ancestor([_ for _ in ids if next(tree.find_elements(name=_, terminal=True), False)])
        while mrca.branch_length == 0:
            mrca = tree.get_path(target=mrca)[-2]
        mrca_ci = format_date_ci(mrca.confidence)
    if mrca != tree.root:
        path = tree.get_path(target=mrca)
        if len(path) == 1:
            parent_ci = format_date_ci(tree.root.branch_length)
        else:
            i = -2
            while i > -len(path) and path[i].branch_length == 0:
                i -= 1
            parent_ci = format_date_ci(tree.root.branch_length if i == -len(path) else path[i].confidence)
    else:
        parent_ci = ''
    print(mrca_ci, parent_ci)
    with open(params.dates, 'w') as f:
        f.write('{}\t{}\t{}\t{}\t{}'.format(params.label, params.subtype, mrca_ci, parent_ci, rate))
