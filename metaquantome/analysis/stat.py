import numpy as np
from scipy import stats as sps
from statsmodels.sandbox.stats import multicomp as mc

from metaquantome.SampleGroups import SampleGroups
from metaquantome.util.constants import P_COLNAME, P_CORR_COLNAME
from metaquantome.util.stat_io import read_expanded_table, write_test


def stat(infile, samps, paired, parametric, ontology, mode, outfile):
    # todo: doc

    # define sample groups
    samp_grps = SampleGroups(samps)

    # read in
    df = read_expanded_table(infile, samp_grps)

    if samp_grps.ngrps != 2:
        ValueError('testing is only available for 2 experimental conditions.')
    # run test
    df_test = test_norm_intensity(df, samp_grps, paired, parametric)
    # write out
    if outfile:
        write_test(df_test, samp_grps=samp_grps, ontology=ontology, mode=mode, outfile=outfile)
    # return
    return df_test


def test_norm_intensity(df, samp_grps, paired, parametric):
    """
    todo: doc
    :param parametric:
    :param df: intensity df, missing values are NaN
    :param samp_grps: is a SampleGroups() object
    :param paired:
    :return:
    """

    grp1_intcols = samp_grps.sample_names[samp_grps.grp_names[0]]
    grp2_intcols = samp_grps.sample_names[samp_grps.grp_names[1]]

    # change any zeros back to NaN
    df.replace(0, np.nan, inplace=True)

    all_intcols = grp1_intcols + grp2_intcols

    test_df = df

    # test, using logged df
    if parametric:
        # test, using logged df
        test_results = test_df.apply(lambda x: sps.stats.ttest_ind(x[grp1_intcols].dropna(),
                                                                   x[grp2_intcols].dropna(),
                                                                   equal_var=paired).pvalue, axis=1)
    else:
        if paired:
            test_results = test_df.apply(lambda x: sps.wilcoxon(x[grp1_intcols].dropna(),
                                                                x[grp2_intcols].dropna()).pvalue,
                                         axis=1)
        else:
            test_results = test_df.apply(lambda x: sps.ranksums(x[grp1_intcols].dropna(),
                                                                x[grp2_intcols].dropna()).pvalue,
                                         axis=1)

    df_means = fold_change(df, samp_grps, log=True)

    df_means[P_COLNAME] = test_results
    df_means[P_CORR_COLNAME] = mc.fdrcorrection0(test_results, method='indep')[1]
    # reset the index to be 'id' - this is mostly for testing, and doesn't affect the output file
    df_means.set_index('id', drop=False, inplace=True)
    return df_means


def fold_change(df, samp_grps, log=False):
    """
    todo: doc
    :param df:
    :param samp_grps:
    :param log:
    :return:
    """
    mean1 = samp_grps.mean_names[0]
    mean2 = samp_grps.mean_names[1]
    if log:
        fc = df[mean1] - df[mean2]
    else:
        fc = np.log2(df[mean1]/df[mean2])
    df[samp_grps.fc_name] = fc
    return df