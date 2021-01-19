import random
import re

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


def get_subtype(name):
    st = re.findall(r'^[^.]+', name)[0]
    if st == 'x':
        return None
    if st[0].isdigit() and st[2] == '_':
        st = 'CRF_{}'.format(st).replace('_cpx', '')
    else:
        st = ', '.join('CRF_{}'.format(_) if _[0].isdigit() else _
                       for _ in re.findall(r'[\d]{2}|A\d|F\d|[ABCDGHJKFUO]', st))
    return st


def get_id(long_id):
    return re.search(r'[^.]+$', long_id)[0]


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_aln', required=True, type=str)
    parser.add_argument('--out_aln', required=True, type=str)
    params = parser.parse_args()

    SeqIO.write((SeqIO.SeqRecord(id=get_id(seq.id), seq=seq.seq, description='')
                 for seq in SeqIO.parse(params.in_aln, 'fasta', alphabet=generic_dna)), params.out_aln, "fasta")


