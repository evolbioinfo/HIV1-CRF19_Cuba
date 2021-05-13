import random

import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

'''
Breakpoints (inclusive, starting from 1) from Casado et al. [2005] with respect to the reference sequence, 
taken from https://www.hiv.lanl.gov/content/sequence/HIV/CRFs/breakpoints.html#CRF19
'''

CRF_19 = 'CRF_19'

subtype2outgroup = {'A1': 'A6', 'D': 'B', 'G': 'A', 'CRF_19': 'C'}


def cut_sequence(recs, breakpoints):
    for rec in recs:
        seq = str(rec.seq)
        cut_seq = ''.join(seq[bp[0]: bp[1]] for bp in breakpoints)
        if len(cut_seq.replace('-', '')) > 100:
            yield SeqIO.SeqRecord(id=rec.id, seq=Seq(cut_seq), description='')


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_aln', required=True, type=str)
    parser.add_argument('--out_aln', required=True, type=str)
    parser.add_argument('--in_data', required=True, type=str)
    parser.add_argument('--crf19_breakpoints', required=True, type=str)
    parser.add_argument('--out_outgroup', required=False, type=str, default=None)
    parser.add_argument('--ingroup', required=True, type=str, choices=('A1', 'G', 'D'))
    params = parser.parse_args()

    df = pd.read_csv(params.in_data, sep='\t', index_col=0, header=0)

    crf_19_breakpoints = {}
    with open(params.crf19_breakpoints, 'r') as f:
        for line in f.readlines():
            line = line.strip('\n').strip()
            if line:
                line = line.split('\t')
                crf_19_breakpoints[line[0]] = [[int(_) for _ in bp.split(' ')] for bp in line[1:]]

    recs = []
    outgroup_recs = []
    outgroup = subtype2outgroup[params.ingroup] if params.out_outgroup else None
    for rec in SeqIO.parse(params.in_aln, 'fasta', alphabet=generic_dna):
        subtype = df.loc[rec.id, 'subtype']
        if subtype in [CRF_19, params.ingroup]:
            recs.append(rec)
        elif df.loc[rec.id, 'subtype_CRF_19_{}'.format(params.ingroup)] == params.ingroup:
            recs.append(rec)
        elif subtype and outgroup == subtype and df.loc[rec.id, 'source'] == 'LA':
            outgroup_recs.append(rec)

    if outgroup_recs:
        outgroup_recs = random.sample(outgroup_recs, 5)

    if params.ingroup != 'CRF_19':
        SeqIO.write(cut_sequence(recs + outgroup_recs, crf_19_breakpoints[params.ingroup]), params.out_aln, "fasta")
    else:
        SeqIO.write(recs + outgroup_recs, params.out_aln, "fasta")

    if outgroup:
        with open(params.out_outgroup, 'w+') as f:
            f.write('\n'.join(rec.id for rec in outgroup_recs))


