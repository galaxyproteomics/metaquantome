import numpy as np
from scipy import stats as sps
from statsmodels.sandbox.stats import multicomp as mc

from metaquantome.classes.SampleGroups import SampleGroups
from metaquantome.util.utils import P_COLNAME, P_CORR_COLNAME
from metaquantome.util.stat_io import read_expanded_table, write_stat


def stat(infile, sinfo, paired, parametric, ontology, mode, outfile):
    """
    Module function that tests differential expression between 2 experimental conditions

    :param infile: path to filtered file
    :param sinfo: path to experimental design file or JSON string
    :param paired: Whether or not the sample should be analyzed as paired samples
    :param parametric: whether or not parameteric tests should be used
    :param ontology: for function, is either 'go', 'ec', or 'cog'
    :param mode: 't', 'f', or 'taxf'
    :param outfile: path to write to
    :return: original dataframe with p value and fold change columns appended
    """

    # define sample groups
    samp_grps = SampleGroups(sinfo)

    # read in
    df = read_expanded_table(infile, samp_grps)

    if samp_grps.ngrps != 2:
        ValueError('testing is only available for 2 experimental conditions.')
    # run test
    df_test = test_norm_intensity(df, samp_grps, paired, parametric)
    # write out
    if outfile:
        write_stat(df_test, samp_grps=samp_grps, ontology=ontology, mode=mode, outfile=outfile)
    # return
    return df_test


def test_norm_intensity(df, samp_grps, paired, parametric):
    """
    run t-tests (or nonparametric tests) on dataframe.
    :param df: intensity df, missing values are NaN
    :param samp_grps: is a SampleGroups() object
    :param paired: whether or not to use a paired test
    :param parametric: Whether or not to use a parametric test
    :return: dataframe with appended pvalue columns
    """

    grp1_intcols = samp_grps.sample_names[samp_grps.grp_names[0]]
    grp2_intcols = samp_grps.sample_names[samp_grps.grp_names[1]]

    # change any zeros back to NaN
    df.replace(0, np.nan, inplace=True)

    # make copy, so df isn't changed by this function
    test_df = df.copy()

    # test, using logged df
    if parametric:
        # don't need to split into paired/unpaired, because both are available in ttest_ind
        test_results = test_df.apply(lambda x: sps.stats.ttest_ind(x[grp1_intcols].dropna(),
                                                                   x[grp2_intcols].dropna(),
                                                                   equal_var=paired).pvalue,
                                     axis=1)
    else:
        if paired:
            # wilcoxon is non-parametric equivalent of paired t-test
            test_results = test_df.apply(lambda x: sps.wilcoxon(x[grp1_intcols].dropna(),
                                                                x[grp2_intcols].dropna()).pvalue,
                                         axis=1)
        else:
            # rank sum test is nonparametric equivalent of unpaired t-test
            test_results = test_df.apply(lambda x: sps.ranksums(x[grp1_intcols].dropna(),
                                                                x[grp2_intcols].dropna()).pvalue,
                                         axis=1)

    # append fold changes to df
    df_means = log2_fold_change(df, samp_grps)

    # p values, uncorrected for multiple comparisons
    df_means[P_COLNAME] = test_results

    # fdr correction
    df_means[P_CORR_COLNAME] = mc.fdrcorrection0(test_results, method='indep')[1]

    # reset the index to be 'id' - this is mostly for testing, and doesn't affect the output file
    df_means.set_index('id', drop=False, inplace=True)
    return df_means


def log2_fold_change(df, samp_grps):
    """
    calculate fold change - fixed as samp_grps.mean_names[0] over samp_grps.mean_names[1],
    where the mean names are sorted alphabetically. The log has already been taken,
    so the L2FC is calculated as mean0 - mean1

    :param df: expanded and/or filtered dataframe
    :param samp_grps: SampleGroups() object
    :return: dataframe with fold change column appended, with name as in samp_grps.fc_name
    """
    mean1 = samp_grps.mean_names[0]
    mean2 = samp_grps.mean_names[1]
    df[samp_grps.fc_name] = df[mean1] - df[mean2]
    return df
