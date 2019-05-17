import os
import pathlib
import numpy as np
import sys
import scipy.stats
import re
import os
import pandas as pd
import seaborn as sns

from VMSfunctions.Controller import *


def get_N(row):
    if 'T10' in row['filename']:
        return 10
    else:
        return row['filename'].split('_')[3]


def experiment_group(row):
    if 'experiment' in row:
        col_to_check = 'experiment'
    else:
        col_to_check = 'filename'

    if 'beer' in row[col_to_check]:
        return 'beer'
    else:
        return 'urine'


def add_group_column(df):
    df['group'] = df.apply(lambda row: experiment_group(row), axis=1)


def get_df(csv_file, min_ms1_intensity, rt_range, mz_range):
    df = pd.read_csv(csv_file)
    intensity_col = 'maxo'
    df = df[(df['rt'] > rt_range[0][0]) & (df['rt'] < rt_range[0][1])]
    df = df[(df['rt'] > mz_range[0][0]) & (df['rt'] < mz_range[0][1])]
    df = df[(df[intensity_col] > min_ms1_intensity)]
    # add log intensity column
    df['log_intensity'] = df.apply(lambda row: np.log(row[intensity_col]), axis=1)
    # add N column
    try:
        df['N'] = df.apply(lambda row: get_N(row), axis=1)
        df[['N']] = df[['N']].astype('int')
    except IndexError:
        pass
    # add group column
    df['group'] = df.apply(lambda row: experiment_group(row), axis=1)
    return df


def make_boxplot(df, x, y, xticklabels, title):
    g = sns.catplot(x=x, y=y, kind='box', data=df)
    g.fig.set_size_inches(10, 3)
    if xticklabels is not None:
        g.set_xticklabels(xticklabels, rotation=90)
    else:
        g.set_xticklabels(rotation=90)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def make_hist(df, col_name, file_name, title):
    gb = df.groupby('filename')
    group_df = gb.get_group(file_name)
    vals = group_df[col_name].values
    print(vals, len(vals))
    _ = plt.hist(vals, bins=100)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def to_chemical(row):
    mz = row['mz'] - PROTON_MASS
    rt = row['rt']
    max_intensity = row['maxo']
    chrom = None
    chem = UnknownChemical(mz, rt, max_intensity, chrom, children=None)
    return chem


def df_to_chemicals(df, filename):
    filtered_df = df.loc[df['filename'] == filename]
    chems = filtered_df.apply(lambda row: to_chemical(row), axis=1).values
    return chems


def find_chem(to_find, min_rts, max_rts, min_mzs, max_mzs, chem_list):
    query_mz = to_find.isotopes[0][0]
    query_rt = to_find.rt
    min_rt_check = min_rts <= query_rt
    max_rt_check = query_rt <= max_rts
    min_mz_check = min_mzs <= query_mz
    max_mz_check = query_mz <= max_mzs
    idx = np.nonzero(min_rt_check & max_rt_check & min_mz_check & max_mz_check)[0]
    matches = chem_list[idx]

    # pick a match
    if len(matches) == 0:
        return None
    elif len(matches) == 1:
        return matches[0]
    else: # multiple matches, take the closest in intensity
        intensity_diffs = [np.abs(chem.max_intensity - to_find.max_intensity) for chem in matches]
        idx = np.argmin(intensity_diffs)
        return matches[idx]


def match(chemical_list_1, chemical_list_2, mz_tol, rt_tol, verbose=False):
    matches = {}
    chem_list = np.array(chemical_list_2)
    min_rts = np.array([chem.rt - rt_tol for chem in chem_list])
    max_rts = np.array([chem.rt + rt_tol for chem in chem_list])
    min_mzs = np.array([chem.isotopes[0][0] * (1-mz_tol/1e6) for chem in chem_list])
    max_mzs = np.array([chem.isotopes[0][0] * (1+mz_tol/1e6) for chem in chem_list])
    for i in range(len(chemical_list_1)):
        to_find = chemical_list_1[i]
        if i % 1000 == 0 and verbose:
            print('%d/%d found %d' % (i, len(chemical_list_1), len(matches)))
        match = find_chem(to_find, min_rts, max_rts, min_mzs, max_mzs, chem_list)
        if match:
            matches[to_find] = match
    return matches
