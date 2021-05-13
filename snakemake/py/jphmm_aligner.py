import re
from collections import defaultdict

from Bio import SeqIO
from Bio.Alphabet import generic_dna
import numpy as np


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_fa', nargs='+', type=str)
    parser.add_argument('--in_msa', nargs='+', type=str)
    parser.add_argument('--out_aln', required=True, type=str)
    parser.add_argument('--aln_len', required=True, type=int)
    params = parser.parse_args()

    id2pos = defaultdict(list)
    for in_msa in params.in_msa:
        with open(in_msa, 'r') as f:
            id = None
            for line in f.readlines():
                line = line.strip('\n').strip()
                if line.startswith('#') or not line:
                    continue
                if line.startswith('>'):
                    id = line[1:].strip()
                    continue
                id2pos[id].extend((int(_) - 1) for _ in line.split(' '))

    with open(params.out_aln, 'w') as f:
        for in_fa in params.in_fa:
            for rec in SeqIO.parse(in_fa, 'fasta', alphabet=generic_dna):
                positions = id2pos[rec.id]
                if positions:
                    aligned_seq = np.array(['-'] * params.aln_len)
                    seq = list(str(rec.seq))
                    insertion_beginning = sum(1 for _ in positions if _ == -1)
                    insertion_end = sum(1 for _ in positions if _ == params.aln_len)
                    aligned_seq[positions[insertion_beginning: (len(positions) - insertion_end)]] \
                        = seq[insertion_beginning: (len(seq) - insertion_end)]
                    f.write('>{}\n{}\n'.format(rec.id, ''.join(aligned_seq)))
