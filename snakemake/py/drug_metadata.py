import os
import re

import pandas as pd
from ete3 import Tree
from sierrapy import SierraClient

QUERY = '''
drugResistance {
    gene { name },
    drugScores {
        drugClass { name },
        drug { displayAbbr, fullName },
        score,
        text,
        partialScores {
            mutations {
                text,
            },
        },
    }
}
'''


def read_tree(nwk):
    tree = None
    for format in range(10):
        try:
            tree = Tree(nwk, format=format)
            break
        except:
            continue
    return tree


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

    parser.add_argument('--input_tab', required=True, type=str)
    parser.add_argument('--output_drug_tab', required=True, type=str)
    parser.add_argument('--output_tab', required=True, type=str)
    parser.add_argument('--timeline', required=True, type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.input_tab, index_col=0, sep='\t')

    data = []
    drms = [str(col) for col in df.columns if col.startswith('RT:') or col.startswith('PR:')]

    drug_bk = params.output_drug_tab + '.backup'
    if os.path.exists(drug_bk):
        drug_df = pd.read_csv(drug_bk, sep='\t', header=0)
        drug_df['mutation'] = drug_df['mutation'].apply(lambda _: str(_).replace('RT_', 'RT:').replace('PR_', 'PR:'))
    else:
        # query Sierra in batches as there is a certain size limit
        for drm in drms:
            for dr in SierraClient().mutations_analysis([drm], QUERY)['drugResistance']:
                gene = dr['gene']['name']
                drug_scores = dr['drugScores']
                for ds in drug_scores:
                    text = ds['text']
                    if 'Resistance' not in text:
                        continue
                    drug_class = ds['drugClass']['name']
                    drug_abbr = ds['drug']['displayAbbr'].replace('/r', '')
                    drug_name = ds['drug']['fullName'].replace('/r', '')
                    score = ds['score']
                    partial_scores = ds['partialScores']
                    for ps in partial_scores:
                        for _ in ps['mutations']:
                            drm = _['text']
                            data.append(['{}:{}'.format(gene, drm), drug_class, drug_name, drug_abbr, score, text])
    
        drug_df = pd.DataFrame(data,
                               columns=['mutation', 'drug class', 'drug full name', 'drug abbreviation', 'score', 'note'])
        drug_df.drop_duplicates(inplace=True)
    
        timeline_df = pd.read_csv(params.timeline, sep='\t', header=0)[['drug abbreviation', 'year']]
        print(timeline_df.head())
        print(drug_df.head())
        drug_df = drug_df.merge(timeline_df, on=['drug abbreviation'], how='inner')
    
        mutations = list(drug_df['mutation'].unique())
        for mut in mutations:
            if not mut.startswith('PR:') and not mut.startswith('RT:') and not mut.startswith('IN:'):
                drug_df = drug_df[drug_df['mutation'] != mut]
                continue
            to_letters = re.sub('[A-Z]\d+', '', mut[mut.find(':') + 1:])
            if len(to_letters) > 1 and to_letters != 'Insertion':
                mut_part = drug_df[drug_df['mutation'] == mut]
                mut = mut.replace(to_letters, '')
                for letter in to_letters:
                    drug_df = drug_df.append(mut_part.copy().assign(**{'mutation': mut + letter}), ignore_index=True)
                drug_df = drug_df[drug_df['mutation'] != mut + to_letters]

    drugs = drug_df['drug abbreviation'].unique()
    for drug in drugs:
        df[drug] = 'sensitive'
        mutations = [_ for _ in drug_df.loc[drug_df['drug abbreviation'] == drug, 'mutation'].unique()
                     if _ in df.columns]
        for mut in mutations:
            df.loc[pd.isna(df[mut]), drug] = None
        for mut in mutations:
            df.loc[df[mut] == 'resistant', drug] = 'resistant'

    df.loc['sensitive', list(drms) + list(drugs)] = 'sensitive'
    df.columns = [_.replace(':', '_') for _ in df.columns]
    drug_df['mutation'] = drug_df['mutation'].apply(lambda _: str(_).replace(':', '_'))
    drug_df.to_csv(params.output_drug_tab, sep='\t', index=False)
    if not os.path.exists(drug_bk):
        drug_df.to_csv(drug_bk, sep='\t', index=False)
    df.to_csv(params.output_tab, sep='\t', index_label='sample_id')
