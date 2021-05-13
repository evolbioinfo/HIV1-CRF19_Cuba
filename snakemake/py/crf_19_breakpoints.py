import random
import re

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

'''
Breakpoints (inclusive, starting from 1) from Casado et al. [2005] with respect to the reference sequence, 
taken from https://www.hiv.lanl.gov/content/sequence/HIV/CRFs/breakpoints.html#CRF19
'''

subtype2breakpoints = {'A1': [(776, 1207), (4240, 4670), (5011, 5290), (6416, 8190), (8511, 8720)],
                       'D': [(1208, 4239), (8721, 9412)],
                       'G': [(4671, 5010), (5291, 6415), (8191, 8510)]}
REFERENCE = 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455'


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_aln', required=True, type=str)
    parser.add_argument('--out_log', required=True, type=str)
    params = parser.parse_args()

    ref = None
    for rec in SeqIO.parse(params.in_aln, 'fasta', alphabet=generic_dna):
        if REFERENCE == rec.id:
            ref = rec
            break

    gappy_seq = str(ref.seq)
    seq = gappy_seq.replace('-', '')

    subtype2gappy_breakpoints = {}
    for subtype, breakpoints in subtype2breakpoints.items():
        shift = 0
        # old breakpoints start from 1 and are inclusive: []
        # new breakpoints start from 0 and are non-inclusive in the end: [)
        gappy_breakpoints = []
        cur_gappy_bp_start = None
        bp_i = 0
        for i in range(0, len(seq)):
            while gappy_seq[i + shift] != seq[i]:
                shift += 1
            if i == breakpoints[bp_i][0] - 1:
                cur_gappy_bp_start = i + shift
            elif i == breakpoints[bp_i][1]:
                gappy_breakpoints.append((cur_gappy_bp_start, i + shift))
                bp_i += 1
                if bp_i == len(breakpoints):
                    break
        subtype2gappy_breakpoints[subtype] = gappy_breakpoints

    with open(params.out_log, 'w') as f:
        for subtype, gappy_breakpoints in subtype2gappy_breakpoints.items():
            f.write('{}\t{}\n'.format(subtype, '\t'.join('{} {}'.format(*bp) for bp in gappy_breakpoints)))


