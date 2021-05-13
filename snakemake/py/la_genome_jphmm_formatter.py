import re
from collections import defaultdict

from Bio import SeqIO
from Bio.Alphabet import generic_dna


def get_subtype(name):
    st = re.findall(r'^[^.]+', name)[0]
    if st == 'x' or st[0].isdigit():
        return []
    return ['CRF_{}'.format(_) if _[0].isdigit() else _ for _ in re.findall(r'[\d]{2}|A\d|F\d|[ABCDGHJKFUO]', st)]


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_aln', required=True, type=str)
    parser.add_argument('--out_aln', required=True, type=str)
    params = parser.parse_args()

    subtype2rec = defaultdict(list)
    for rec in SeqIO.parse(params.in_aln, 'fasta', alphabet=generic_dna):
        sts = get_subtype(rec.id)
        if len(sts) == 1 and sts[0] != 'A':
            subtype2rec[sts[0]].append(rec)
            print(sts)
    with open(params.out_aln, 'w') as f:
        for subtype in sorted(subtype2rec.keys()):
            f.write('>>{}\n'.format(subtype))
            for rec in subtype2rec[subtype]:
                f.write('>{}\n{}\n'.format(rec.id, rec.seq))
