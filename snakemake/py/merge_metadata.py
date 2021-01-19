import pandas as pd
from hdx.location.country import Country


def add_loc_info(df):
    ccs = df['country_code'].unique()

    iso22info = {_: Country.get_country_info_from_iso2(_) for _ in ccs if not pd.isna(_)}

    df['country'] = df['country_code'].apply(lambda _: iso22info[_]['Country or Area'] if _ in iso22info else None)
    df['region'] = df['country_code'].apply(lambda _: iso22info[_]['Region Name'] if _ in iso22info else None)
    df['subregion'] = df['country_code'].apply(lambda _: iso22info[_]['Sub-region Name'] if _ in iso22info else None)
    df['intregion'] = df['country_code'].apply(lambda _: iso22info[_]['Intermediate Region Name'] if _ in iso22info else None)

    asia = df['region'] == 'Asia'
    df.loc[asia, 'intregion'] = df.loc[asia, 'subregion']
    df.loc[asia, 'subregion'] = 'Asia'

    europe = df['region'] == 'Europe'
    df.loc[europe, 'intregion'] = df.loc[europe, 'subregion']
    df.loc[europe, 'subregion'] = 'Europe'

    russia = df['country'] == 'Russian Federation'
    df.loc[russia, 'intregion'] = 'Russia'
    df.loc[russia, 'subregion'] = 'Russia'

    la_and_carrib = df['subregion'] == 'Latin America and the Caribbean'
    df.loc[la_and_carrib, 'subregion'] = df.loc[la_and_carrib, 'intregion']

    cuba = df['country'] == 'Cuba'
    df.loc[cuba, 'intregion'] = 'Cuba'
    df.loc[cuba, 'subregion'] = 'Cuba'
    df.loc[cuba, 'region'] = 'Cuba'

    northen_am = df['subregion'] == 'Northern America'
    df.loc[northen_am, 'intregion'] = 'Northern America'

    aus_nz = df['subregion'] == 'Australia and New Zealand'
    df.loc[aus_nz, 'intregion'] = df.loc[aus_nz, 'country']

    return df


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--data', nargs='+', type=str)
    parser.add_argument('--output_data', required=True, type=str)
    parser.add_argument('--subtype_data', nargs='*', type=str, default=[])
    params = parser.parse_args()

    df = pd.concat([pd.read_csv(tab, sep='\t', index_col=0, header=0) for tab in params.data], join='outer', sort=False)
    df = df.loc[~df.index.duplicated(keep='first')]

    for file in params.subtype_data:
        rec_df = pd.read_table(file, skiprows=9, header=None, index_col=0, names=['params', 'jpHMM subtype'])
        rec_df.index = rec_df.index.str.replace('^>', '')
        df.loc[rec_df.index, 'jpHMM_subtype'] = rec_df['jpHMM subtype']

    # df.loc[~pd.isna(df['country_code']) & (df['country_code'] != 'CU'), 'region_of_residence'] = 'abroad'
    df.loc[~pd.isna(df['country_code']) & (df['country_code'] != 'CU'), 'province_of_diagnostics'] = 'abroad'

    add_loc_info(df)

    df.to_csv(params.output_data, sep='\t', header=True, index=True, index_label='id')
