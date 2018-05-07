import pandas as pd
import numpy as np
import scipy.stats as sps
import statsmodels.sandbox.stats.multicomp as mc

def define_intensity_columns(grp1_intcols, grp2_intcols):
    if not isinstance(grp1_intcols, list):
        grp1 = [grp1_intcols]
    else:
        grp1 = grp1_intcols

    if not isinstance(grp2_intcols, list):
        grp2 = [grp2_intcols]
    else:
        grp2 = grp2_intcols

    if grp2_intcols:
        all_intcols = grp1 + grp2
    else:
        all_intcols = grp1

    dict_numeric_cols = {x: np.float64 for x in all_intcols}

    return all_intcols, dict_numeric_cols


def read_data_table(file, dict_numeric_cols, all_intcols, pep_colname="peptide"):
    # read in data
    df = pd.read_table(file, sep="\t", index_col=pep_colname, dtype=dict_numeric_cols,
                       na_values = ["", "0", "NA", "NaN"])

    # drop columns where all are NA
    df.dropna(axis=1, how="all", inplace=True)

    # change missing intensities to 0
    values = {x: 0 for x in all_intcols}
    df.fillna(values, inplace=True)

    # then change missing taxa/function to 'unknown'
    df.fillna('unknown', inplace=True)

    return df


def filter_min_observed(df, grp1_intcols, grp2_intcols, threshold):
    if threshold == 0:
        return df
    # filter to a minimum number of observed intensities per sample
    numNotNA1 = (df[grp1_intcols] > 0).apply(sum, axis=1) >= threshold
    numNotNA2 = (df[grp2_intcols] > 0).apply(sum, axis=1) >= threshold
    keep = numNotNA1 & numNotNA2
    filtered_df = df.loc[keep].copy()
    return filtered_df


def test_norm_intensity(df, grp1_intcols, grp2_intcols, paired):

    # add mean, sd; create dataframe
    test_results = df.apply(lambda x: sps.stats.ttest_ind(x[grp1_intcols],
                                                                   x[grp2_intcols],
                                                                   equal_var=paired).pvalue, axis=1)
    df['log2ratio_2over1'] = \
        df.apply(
            lambda x: np.log2(np.mean(x[grp2_intcols]) / np.mean(x[grp1_intcols])),
        axis=1)
    df['p'] = test_results
    df['corrected_p'] = mc.fdrcorrection0(test_results, method='indep')[1]
    df['id'] = df.index
    return df