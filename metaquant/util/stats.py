import numpy as np
import scipy.stats as sps
import statsmodels.sandbox.stats.multicomp as mc

P_COLNAME = 'p'
P_CORR_COLNAME = 'corrected_p'


def group_and_sum_by_rank(df, rank, all_intcols, norm_to_rank=False):
    # sum intensities in each rank
    summed_abund = df.groupby(by=rank)[all_intcols].sum(axis=0)

    # normalize to each sample - not currently done
    if norm_to_rank:
        return_df = summed_abund / summed_abund.sum(axis=0)
    else:
        return_df = summed_abund

    return_df['rank'] = rank
    return_df['id'] = return_df.index
    return return_df


def test_norm_intensity(df, samp_grps, paired, parametric):
    """

    :param parametric:
    :param df:
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

    return df_means


def calc_means(df, samp_grps):

    for i in range(samp_grps.ngrps):
        grp_name = samp_grps.grp_names[i]
        samples_in_grp = samp_grps.sample_names[grp_name]
        if len(samples_in_grp) > 1:
            sample = df[samples_in_grp]
            means = np.mean(sample, axis = 1)
            means[means == 0] = np.nan
            logs = np.log2(means)
            df[samp_grps.mean_names[i]] = logs

        else:
            # just log transform the single sample
            df[samp_grps.mean_names[i]] = np.log2(df[samples_in_grp])

    return df


def fold_change(df, samp_grps, log=False):
    mean1 = samp_grps.mean_names[0]
    mean2 = samp_grps.mean_names[1]
    if log:
        fc = df[mean1] - df[mean2]
    else:
        fc = np.log2(df[mean1]/df[mean2])
    df[samp_grps.fc_name] = fc
    return df
