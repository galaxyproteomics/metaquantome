from metaquant.AnnotationHierarchy import AnnotationHierarchy
from metaquant.SampleAnnotations import SampleAnnotations
from metaquant.util import stats
from metaquant.databases.NCBITaxonomyDb import NCBITaxonomyDb
import numpy as np


def common_hierarchical_analysis(db, df, annot_colname, samp_grps,
                                 min_peptides, min_children_non_leaf,
                                 test, threshold, paired, parametric):
    samp_annot = SampleAnnotations(db)
    # make a hierarchy for each sample
    samp_annot.add_samples_from_df(df, annot_colname, samp_grps, min_peptides, min_children_non_leaf)
    intensity_all_ranks = samp_annot.to_dataframe()

    # change zeros back to NaN for stats, so they are ignored
    results = common_stats(intensity_all_ranks, samp_grps, test, threshold, paired, parametric)
    results['id'] = results.index
    return results


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