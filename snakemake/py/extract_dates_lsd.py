import pandas as pd
from ete3 import Tree

MIN_DATE = 1900
MAX_DATE = 2019


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


def get_tip_name(n):
    return n.name if n.is_leaf() else next(_.name for _ in n)


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--data', required=True, type=str)
    parser.add_argument('--dates', required=True, type=str)
    parser.add_argument('--tree', required=True, type=str)
    parser.add_argument('--date_col', required=True, type=str)
    parser.add_argument('--diag_col', required=True, type=str)
    params = parser.parse_args()

    # we add additional constraints that ensure that if A and B have diagnostics dates that are older
    # than their corresponding sampling dates,
    # then their MRCA's date cannot be after the last diagnosis date (i.e. they should diverge before being diagnosed)
    tree = Tree(params.tree)

    date_df = pd.read_csv(params.data, index_col=0, sep='\t')
    date_df.index = date_df.index.map(str)
    date_df = date_df.loc[[_.name for _ in tree], :]

    date_df[params.date_col] = pd.to_datetime(
        date_df[params.date_col].astype(str).str.replace('.0', '', regex=False),
        infer_datetime_format=True, errors='coerce').map(_get_date)
    if params.diag_col in date_df.columns:
        date_df[params.diag_col] = pd.to_datetime(
            date_df[params.diag_col].astype(str).str.replace('.0', '', regex=False),
            infer_datetime_format=True, errors='coerce').map(_get_date)

        for tip in tree:
            date = date_df.loc[tip.name, params.date_col]
            diag_date = date_df.loc[tip.name, params.diag_col]
            if not pd.isna(date):
                date_df.loc[tip.name, 'lsd_date'] = date
            elif not pd.isna(diag_date):
                date_df.loc[tip.name, 'lsd_date'] = 'b({},{})'.format(diag_date, MAX_DATE)
            else:
                date_df.loc[tip.name, 'lsd_date'] = 'b({},{})'.format(MIN_DATE, MAX_DATE)
            tip.add_feature('diag_date', diag_date)

        for n in tree.traverse('postorder'):
            if not n.is_leaf():
                diag_dates = sorted([getattr(_, 'diag_date') for _ in n.children
                                     if not pd.isna(getattr(_, 'diag_date'))])
                n.add_feature('diag_date', None if not diag_dates else diag_dates[0])
                if len(diag_dates) > 1:
                    date_df.loc['mrca({})'.format(','.join(get_tip_name(c) for c in n.children)), 'lsd_date'] \
                        = 'b({},{})'.format(MIN_DATE, diag_dates[1] - 1 / 2 / 365)
    else:
        for tip in tree:
            date = date_df.loc[tip.name, params.date_col]
            if not pd.isna(date):
                date_df.loc[tip.name, 'lsd_date'] = date

    date_df = date_df[~pd.isna(date_df['lsd_date'])]
    with open(params.dates, 'w+') as f:
        f.write('%d\n' % date_df.shape[0])
    date_df['lsd_date'].to_csv(params.dates, sep='\t', header=False, mode='a')



