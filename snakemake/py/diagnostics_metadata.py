import pandas as pd
import numpy as np


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


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_data', required=True, type=str)
    parser.add_argument('--out_data', required=True, type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.in_data, index_col=0, header=0, sep='\t')
    df = df[df['subtype'].isin(['CRF_19', 'D', 'A1', 'G'])
            | (df['subtype_CRF_19_D'] == 'D')
            | (df['subtype_CRF_19_A1'] == 'A1')
            | (df['subtype_CRF_19_G'] == 'G')]
    df['subtype'] = np.where(df['source'] == 'LA', df['subtype'],
                             np.where(df['subtype_preannotated'] == df['subtype'], df['subtype'], None))
    df['location'] = np.where(df['subregion'] == 'Sub-Saharan Africa', df['country'], df['subregion'])

    if 'sexuality' in df.columns:
        df.loc[df['sexuality'] == 'Bisexual', 'sexuality'] = None

    if 'diagnostics_date' in df.columns:
        diagnostics_columns = [_ for _ in ['reason_for_detection', 'province_of_diagnostics', 'AIDS_at_diagnostics']
                               if _ in df.columns]

        df['sample_date'] = pd.to_datetime(df['sample_date'].astype(str).str.replace('.0', '', regex=False),
                                           infer_datetime_format=True, errors='coerce').map(_get_date)
        df['diagnostics_date'] = pd.to_datetime(df['diagnostics_date'].astype(str).str.replace('.0', '', regex=False),
                                                infer_datetime_format=True, errors='coerce').map(_get_date)

        ids = list(df[~pd.isna(df['sample_date']) & ~pd.isna(df['diagnostics_date'])
                      & (df['sample_date'] > df['diagnostics_date'])].index)

        for id in ids:
            df.loc[id + '_diagnostics', diagnostics_columns] = df.loc[id, diagnostics_columns]
            df.loc[id + '_diagnostics', ['country_code', 'country', 'intregion', 'subregion', 'subtype', 'source', 'location']] \
                = df.loc[id, ['country_code', 'country', 'intregion', 'subregion', 'subtype', 'source', 'location']]
            if 'sexuality' in df.columns:
                df.loc[id + '_diagnostics', 'sexuality'] = df.loc[id, 'sexuality']

        df.loc[ids, diagnostics_columns] = [None] * len(diagnostics_columns)

    df.to_csv(params.out_data, sep='\t', header=True, index=True, index_label='id')
