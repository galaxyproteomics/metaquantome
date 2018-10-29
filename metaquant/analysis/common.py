from metaquant.SampleAnnotations import SampleAnnotations
from metaquant.util import stats
import numpy as np


def common_hierarchical_analysis(db, df, annot_colname, samp_grps,
                                 min_peptides, min_children_non_leaf, threshold):
    samp_annot = SampleAnnotations(db)
    # make a hierarchy for each sample
    samp_annot.add_samples_from_df(df, annot_colname, samp_grps, min_peptides, min_children_non_leaf)
    intensity_all_ranks = samp_annot.to_dataframe()

    # filter
    int_all_ranks_filt = stats.filter_min_observed(intensity_all_ranks, threshold, samp_grps)
    int_all_ranks_filt['id'] = intensity_all_ranks.index

    # calculate means
    int_w_means = stats.calc_means(int_all_ranks_filt, samp_grps)

    # clean and log transform
    # replace nan with zero, so that np.log2 returns nan
    int_w_means[int_w_means == 0] = np.nan
    # take log of intensities for return
    int_w_means[samp_grps.all_intcols] = np.log2(int_w_means[samp_grps.all_intcols])

    return int_w_means

    # old
    # change zeros back to NaN for stats, so they are ignored
    # results = common_stats(intensity_all_ranks, samp_grps, test, threshold, paired, parametric)
    # results['id'] = results.index
    # return results

# todo: move this to stats module
def common_stats(intensity_df, samp_grps, test, threshold, paired, parametric):
    # test
    if test and samp_grps.ngrps == 2:
        results = stats.test_norm_intensity(intensity_df, samp_grps, threshold, paired, parametric)
    else:
        results = stats.calc_means(intensity_df, samp_grps)
    # replace nan with zero, so that np.log2 returns nan
    results[results == 0] = np.nan
    # take log of intensities for return
    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])
    return results
