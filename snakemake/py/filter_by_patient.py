import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna


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
        return pd.np.inf


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--data', required=True, type=str)
    parser.add_argument('--in_aln', required=True, type=str)
    parser.add_argument('--out_aln', required=True, type=str)
    parser.add_argument('--date_col', required=True, type=str)
    parser.add_argument('--patient_col', required=True, type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.data, index_col=0, sep='\t')
    if params.patient_col in df.columns:
        df = df[[params.date_col, params.patient_col]]
        df[params.date_col] = pd.to_datetime(df[params.date_col].astype(str).str.replace('.0', '', regex=False),
                                             infer_datetime_format=True, errors='coerce').map(_get_date)
    df.index = df.index.map(str)

    def filter_recs():
        patients = set()
        for _ in sorted(list(SeqIO.parse(params.in_aln, "fasta", alphabet=generic_dna)),
                        key=lambda _: df.loc[_.id, params.date_col] if _.id in df.index else pd.np.inf):
            if params.patient_col in df.columns and _.id in df.index:
                patient = df.loc[_.id, params.patient_col]
                if pd.isna(patient):
                    yield _
                elif patient not in patients:
                    patients.add(patient)
                    yield _
            else:
                yield _

    count = SeqIO.write(filter_recs(), params.out_aln, "fasta")



