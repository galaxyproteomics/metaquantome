import itertools
import numpy as np


def separation_of_n_clusts(df, samp_grps):
    # number of experimental conditions
    ngrps = samp_grps.ngrps

    # get the center of each cluster in feature space
    centroids = centroid(df, samp_grps)

    # average distance between centroids
    avg_dist_bet_cent = avg_dist_between_group(centroids)

    # within variance
    tot_within = total_sse_within_group(df, centroids, samp_grps)

    return avg_dist_bet_cent / tot_within


def avg_dist_between_group(centroids):
    # possible combinations of conditions
    combs = list(itertools.combinations(centroids, 2))
    # between variance
    # iterate over combinations
    ndists = len(combs)
    dists = [0] * ndists
    for i in range(ndists):
        ith_comb = combs[i]
        dists[i] = np.square((ith_comb[0] - ith_comb[1])).values.sum()
    avg_dist = np.mean(dists)
    return avg_dist


def centroid(df, samp_grps):
    """
    calculate the mean of each condition in feature space
    :param df: DataFrame
    :param samp_grps: SampleGroups object
    :return: list of Series in the order of samp_grps.grp_names
    """
    # this takes the row means for each sample group
    means = [df[samp_grps.sample_names[elem]].mean(axis=1) for elem in samp_grps.grp_names]
    return means


def total_sse_within_group(df, centroids, samp_grps):
    """
    Using a list of centroids, calculate the average variance within each group
    :param df: DataFrame with intensities for each term
    :param means:
    :param samp_grps:
    :return:
    """
    # ngroups
    ngrps = samp_grps.ngrps
    # within variance
    within = [0] * ngrps
    for i in range(ngrps):
        grp = samp_grps.grp_names[i]
        within[i] = np.square(df[samp_grps.sample_names[grp]].subtract(centroids[i], axis='index')).values.sum()
    tot_within = np.sum(within)
    return tot_within

