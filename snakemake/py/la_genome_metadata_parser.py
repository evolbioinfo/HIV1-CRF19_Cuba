import re

import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from hdx.location.country import Country

subtypes = {'A1', 'D', 'G', 'CRF_19'}


def get_date(name):
    # rm subtype
    name = re.sub(r'^[^.]+\.', '', name)
    # rm country code
    name = re.sub(r'^[^.]+\.', '', name)
    date = re.findall(r'^[\d]{2}|x', name)[0]
    if date == 'x':
        return None
    if date.startswith('0') or date.startswith('1'):
        return '20' + date
    return '19' + date


def get_cc(name):
    # rm subtype
    name = re.sub(r'^[^.]+\.', '', name)
    loc = re.findall(r'^[\w]{2}|x', name)[0]
    if loc == 'x':
        return None
    return loc


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


def add_loc_info(df):
    ccs = df['country_code'].unique()

    iso22info = {_: Country.get_country_info_from_iso2(_) for _ in ccs if _}

    df['country'] = df['country_code'].apply(lambda _: iso22info[_]['#country+name+preferred'] if _ in iso22info else None)
    df['region'] = df['country_code'].apply(lambda _: iso22info[_]['#region+main+name+preferred'] if _ in iso22info else None)
    df['sub-region'] = df['country_code'].apply(lambda _: iso22info[_]['#region+name+preferred+sub'] if _ in iso22info else None)
    df['int-region'] = df['country_code'].apply(lambda _: iso22info[_]['#region+intermediate+name+preferred'] if _ in iso22info else None)

    asia = df['region'] == 'Asia'
    df.loc[asia, 'int-region'] = df.loc[asia, 'sub-region']
    df.loc[asia, 'sub-region'] = 'Asia'

    europe = df['region'] == 'Europe'
    df.loc[europe, 'int-region'] = df.loc[europe, 'sub-region']
    df.loc[europe, 'sub-region'] = 'Europe'

    russia = df['country'] == 'Russian Federation'
    df.loc[russia, 'int-region'] = 'Russia'
    df.loc[russia, 'sub-region'] = 'Russia'

    la_and_carrib = df['sub-region'] == 'Latin America and the Caribbean'
    df.loc[la_and_carrib, 'sub-region'] = df.loc[la_and_carrib, 'int-region']

    cuba = df['country'] == 'Cuba'
    df.loc[cuba, 'int-region'] = 'Cuba'
    df.loc[cuba, 'sub-region'] = 'Cuba'

    northen_am = df['sub-region'] == 'Northern America'
    df.loc[northen_am, 'int-region'] = 'Northern America'

    aus_nz = df['sub-region'] == 'Australia and New Zealand'
    df.loc[aus_nz, 'int-region'] = df.loc[aus_nz, 'country']

    return df


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_aln', required=True, type=str)
    parser.add_argument('--out_data', required=True, type=str)
    parser.add_argument('--no_loc_info', action='store_true', default=False)
    parser.add_argument('--any_subtype', action='store_true', default=False)
    params = parser.parse_args()

    index, data = [], []
    for rec in SeqIO.parse(params.in_aln, 'fasta', alphabet=generic_dna):
        subtype = get_subtype(rec.id)
        if subtype in subtypes or params.any_subtype:
            index.append(get_id(rec.id))
            data.append([get_date(rec.id), get_cc(rec.id), get_subtype(rec.id)])

    df = pd.DataFrame(data=data, index=index, columns=['sample_date', 'country_code', 'subtype'])
    df['source'] = 'LA'

    if not params.no_loc_info:
        add_loc_info(df)

    df.to_csv(params.out_data, sep='\t', header=True, index=True, index_label='accession')


