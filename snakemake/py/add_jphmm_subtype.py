from collections import defaultdict

import pandas as pd

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_data', nargs='+', type=str)
    parser.add_argument('--in_msa', nargs='+', type=str)
    parser.add_argument('--in_rec', nargs='+', type=str)
    parser.add_argument('--out_data', required=True, type=str)
    parser.add_argument('--crf19_breakpoints', required=False, type=str, default=None)
    parser.add_argument('--slack', required=False, type=float, default=0.4)
    params = parser.parse_args()

    id2pos = defaultdict(list)
    for in_msa in params.in_msa:
        with open(in_msa, 'r') as f:
            id = None
            for line in f:
                line = line.strip('\n').strip()
                if line.startswith('#') or not line:
                    continue
                if line.startswith('>'):
                    id = line[1:].strip()
                    continue
                id2pos[id].extend((int(_) - 1) for _ in line.split(' '))

    df = pd.concat(pd.read_csv(in_data, sep='\t', header=0, index_col=0) for in_data in params.in_data)
    df.index = df.index.map(str)
    df['subtype_preannotated'] = df['subtype']

    id2rec = defaultdict(lambda: defaultdict(list))
    for in_rec in params.in_rec:
        with open(in_rec, 'r') as f:
            id = None
            for line in f:
                line = line.strip('\n').strip()
                if line.startswith('#') or not line:
                    continue
                if line.startswith('>'):
                    id = line[1:].strip().split(' ')[0]
                    continue
                start, end, subtype = line.split('\t')
                id2rec[id][subtype].append((int(start), int(end)))

    if params.crf19_breakpoints:
        crf_19_st2breakpoints = {}
        with open(params.crf19_breakpoints, 'r') as f:
            for line in f:
                line = line.strip('\n').strip()
                if line:
                    line = line.split('\t')
                    crf_19_st2breakpoints[line[0]] = [[int(_) for _ in bp.split(' ')] for bp in line[1:]]

        crf19_subtypes = set(crf_19_st2breakpoints.keys())

        for id, subtype2pos in id2rec.items():
            subtypes = sorted(subtype2pos.keys())
            subtypes = [_ for _ in subtypes if 'Insertion' not in _]
            df.loc[id, 'subtype_jpHMM'] = ','.join(subtypes)

            # check what we have for each CRF_19 breakpoint
            crf_19_subtype2subtypes = defaultdict(set)
            for crf_19_subtype, crf_19_breakpoints in crf_19_st2breakpoints.items():
                for (bp_start, bp_stop) in crf_19_breakpoints:
                    for subtype, breakpoints in subtype2pos.items():
                        for (start, stop) in breakpoints:
                            start = id2pos[id][start - 1]
                            stop = id2pos[id][stop - 1] + 1
                            if start <= bp_stop and stop >= bp_start \
                                    and (min(bp_stop, stop) - max(bp_start, start)) / (stop - start) > params.slack:
                                crf_19_subtype2subtypes[crf_19_subtype].add(subtype)

            within_breakpoints = True
            for crf_19_st in crf19_subtypes:
                sts = crf_19_subtype2subtypes[crf_19_st]
                if sts - {crf_19_st}:
                    within_breakpoints = False
                df.loc[id, 'subtype_CRF_19_{}'.format(crf_19_st)] = ','.join(sorted(sts))

            df.loc[id, 'CRF_19_compatible'] = within_breakpoints
            df.loc[id, 'subtype'] = 'CRF_19' if within_breakpoints and ((df.loc[id, 'subtype_preannotated'] == 'CRF_19') or (len(subtypes) > 1)) \
                else df.loc[id, 'subtype_jpHMM']

    df.columns = [c.replace(' ', '_') for c in df.columns]
    df[~pd.isna(df['subtype_jpHMM'])].to_csv(params.out_data, sep='\t', index_label='id')
