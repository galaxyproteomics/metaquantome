import metaquant.analysis.filter
from metaquant.SampleAnnotations import SampleAnnotations
from metaquant.util import stats
import numpy as np


def common_hierarchical_analysis(db, df, annot_colname, samp_grps, hierarchical=True):

    # import pdb; pdb.set_trace()
    if hierarchical:
        samp_annot = SampleAnnotations(db)
        # make a hierarchy for each sample
        samp_annot.add_samples_from_df(df, annot_colname, samp_grps)
        intensity_all_ranks = samp_annot.to_dataframe()
    else:
        intensity_all_ranks = df

    # calculate means
    int_w_means = stats.calc_means(intensity_all_ranks, samp_grps)

    # clean and log transform
    # replace nan with zero, so that np.log2 returns nan
    int_w_means[int_w_means == 0] = np.nan

    # take log of intensities for return
    int_w_means[samp_grps.all_intcols] = np.log2(int_w_means[samp_grps.all_intcols])

    int_w_means['id'] = int_w_means.index

    return int_w_means


# todo: move this to stats module
def common_stats(intensity_df, samp_grps, test, paired, parametric):
    # test
    if test and samp_grps.ngrps == 2:
        results = stats.test_norm_intensity(intensity_df, samp_grps, paired, parametric)
    else:
        results = stats.calc_means(intensity_df, samp_grps)
    # replace nan with zero, so that np.log2 returns nan
    results[results == 0] = np.nan
    # take log of intensities for return
    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])
    return results
