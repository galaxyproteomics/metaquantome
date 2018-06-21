import numpy as np
import scipy.stats as sps
import statsmodels.sandbox.stats.multicomp as mc


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

    # change any zeros back to NaN
    df_filt.replace(0, np.nan, inplace=True)

    all_intcols = grp1_intcols + grp2_intcols

    test_df = np.log2(df_filt[all_intcols])

    # test, using logged df
    test_results = test_df.apply(lambda x: sps.stats.ttest_ind(x[grp1_intcols].dropna(),
                                                                   x[grp2_intcols].dropna(),
                                                                   equal_var=paired).pvalue, axis=1)

    df_means = fold_change(calc_means(df_filt, samp_grps), samp_grps, log=True)

    df_means['p'] = test_results
    df_means['corrected_p'] = mc.fdrcorrection0(test_results, method='indep')[1]
    df_means['id'] = df_filt.index

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
            df[samp_grps.mean_names] = np.log2(df[samples_in_grp])

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
