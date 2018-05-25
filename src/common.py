import pandas as pd
import numpy as np
import scipy.stats as sps
import statsmodels.sandbox.stats.multicomp as mc
import itertools


class SampleGroups:
    """
    flatten sample names list
    :param sample_names: dictionary of lists, where the keys are the group names and
    the values are lists of column names within that group
    :return:
    """

    def __init__(self, sample_names):

        # top level dictionary
        self.sample_names = sample_names

        # flatten sample column names
        self.all_intcols = list(itertools.chain(*list(sample_names.values())))

        # for defining pandas column data types on read in
        self.dict_numeric_cols = {x: np.float64 for x in self.all_intcols}

        # number of conditions
        self.ngrps = len(sample_names)

        # name of experimental groups
        # sort alphabetically, so it's deterministic
        self.grp_names = sorted(list(sample_names.keys()))

        # when calculating means, column names for means
        self.mean_names = [grp + "_mean" for grp in self.grp_names]



def read_data_table(file, samp_grps, pep_colname):
    # read in data
    df = pd.read_table(file, sep="\t", index_col=pep_colname, dtype=samp_grps.dict_numeric_cols,
                       na_values = ["", "0", "NA", "NaN"], low_memory=False)

    # drop columns where all are NA
    df.dropna(axis=1, how="all", inplace=True)

    # change missing intensities to 0
    values = {x: 0 for x in samp_grps.all_intcols}
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


def test_norm_intensity(df, samp_grps, threshold, paired, log=False):
    """

    :param df:
    :param samp_grps: is a SampleGroups() object
    :param threshold:
    :param paired:
    :param log:
    :return:
    """
    grp1_intcols = samp_grps.sample_names[samp_grps.grp_names[0]]
    grp2_intcols = samp_grps.sample_names[samp_grps.grp_names[1]]

    df_filt = filter_min_observed(df, grp1_intcols, grp2_intcols, threshold)

    # change zeros back to NaN
    df_filt.replace(0, np.nan, inplace=True)

    # add mean, sd; create dataframe
    test_results = df_filt.apply(lambda x: sps.stats.ttest_ind(x[grp1_intcols].dropna(),
                                                                   x[grp2_intcols].dropna(),
                                                                   equal_var=paired).pvalue, axis=1)

    df_means = fold_change(calc_means(df_filt, samp_grps), samp_grps, log)

    df_means['p'] = test_results
    df_means['corrected_p'] = mc.fdrcorrection0(test_results, method='indep')[1]
    df_means['id'] = df_filt.index

    return df_means


def calc_means(df, samp_grps):

    for i in range(samp_grps.ngrps):
        grp_name = samp_grps.grp_names[i]
        samples_in_grp = samp_grps.sample_names[grp_name]

        if len(samples_in_grp) > 1:
            df[samp_grps.mean_names[i]] = df.apply(lambda x: np.mean(x[samples_in_grp].dropna()), axis=1)
        else:
            df[samp_grps.mean_names[i]] = df[samples_in_grp]

    return df


def fold_change(df, samp_grps, log=False):
    mean1 = samp_grps.mean_names[0]
    mean2 = samp_grps.mean_names[1]

    grp1 = samp_grps.grp_names[0]
    grp2 = samp_grps.grp_names[1]

    if log:
        fc = df[mean1] - df[mean2]
    else:
        fc = np.log2(df[mean1]/df[mean2])

    fc_name = 'log2fc_' + grp1 + '_over_' + grp2

    df[fc_name] = fc
    return df
