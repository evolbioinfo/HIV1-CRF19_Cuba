import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna


def subtype_combiner(col, df1, df2, df, common_index):
    df1[col].fillna('', inplace=True)
    df2[col].fillna('', inplace=True)
    df.loc[common_index, col] = \
        pd.np.where(df1.loc[common_index, col] != df2.loc[common_index, col],
                    (df1.loc[common_index, col] + ',' + df2.loc[common_index, col])
                    .apply(lambda _: ','.join(sorted(set(_.split(','))))),
                    df1.loc[common_index, col])


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_data', nargs='+', type=str)
    parser.add_argument('--out_data', required=True, type=str)
    parser.add_argument('--in_aln', nargs='+', type=str)
    parser.add_argument('--out_aln', required=True, type=str)
    params = parser.parse_args()

    df = None
    for tab in params.in_data:
        cur_df = pd.read_csv(tab, sep='\t', index_col=0, header=0)
        if df is None:
            df = cur_df
        else:
            common_cols = df.columns & cur_df.columns
            common_index = df.index & cur_df.index
            for col in common_cols:
                df.loc[common_index, col] = pd.np.where(pd.isna(df.loc[common_index, col]),
                                                        cur_df.loc[common_index, col], df.loc[common_index, col])
                cur_df.loc[common_index, col] = pd.np.where(pd.isna(cur_df.loc[common_index, col]),
                                                            df.loc[common_index, col], cur_df.loc[common_index, col])
            crf_19_compatibility_mask = df.loc[common_index, 'CRF_19_compatible'] & cur_df.loc[common_index, 'CRF_19_compatible']
            df.loc[common_index, 'subtype'] = pd.np.where(crf_19_compatibility_mask, df.loc[common_index, 'subtype'],
                                                          df.loc[common_index, 'subtype_jpHMM'])
            cur_df.loc[common_index, 'subtype'] = pd.np.where(crf_19_compatibility_mask, cur_df.loc[common_index, 'subtype'],
                                                              cur_df.loc[common_index, 'subtype_jpHMM'])

            merged_df = df.merge(cur_df, how='outer', left_index=True, right_index=True, on=list(common_cols))
            merged_df.loc[common_index, 'CRF_19_compatible'] = crf_19_compatibility_mask
            subtype_combiner('subtype_preannotated', df, cur_df, merged_df, common_index)
            subtype_combiner('subtype_jpHMM', df, cur_df, merged_df, common_index)
            subtype_combiner('subtype', df, cur_df, merged_df, common_index)
            subtype_combiner('subtype_CRF_19_D', df, cur_df, merged_df, common_index)
            subtype_combiner('subtype_CRF_19_A1', df, cur_df, merged_df, common_index)
            subtype_combiner('subtype_CRF_19_G', df, cur_df, merged_df, common_index)

            merged_df.loc[common_index, 'subtype'] = pd.np.where(merged_df.loc[common_index, 'CRF_19_compatible'],
                                                                 'CRF_19', merged_df.loc[common_index, 'subtype'])
            df = merged_df

    df['source'] = 'CU'
    df.to_csv(params.out_data, sep='\t', header=True, index=True, index_label='id')

    id2seq = {}
    for fa in params.in_aln:
        for rec in SeqIO.parse(fa, 'fasta', alphabet=generic_dna):
            if rec.id not in id2seq:
                id2seq[rec.id] = str(rec.seq)
            else:
                seq = pd.np.array(list(id2seq[rec.id]))
                seq[seq == '-'] = pd.np.array(list(str(rec.seq)))[seq == '-']
                id2seq[rec.id] = ''.join(seq)

    with open(params.out_aln, 'w') as f:
        for id, seq in id2seq.items():
            f.write('>{}\n{}\n'.format(id, seq))
