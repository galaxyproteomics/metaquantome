from metaquant.AnnotationHierarchy import AnnotationHierarchy
from metaquant.util import stats
from metaquant.databases.NCBITaxonomyDb import NCBITaxonomyDb
import numpy as np


def common_hierarchical_analysis(db, df, annot_colname, samp_grps, min_peptides, min_children_non_leaf, test,
                                 threshold, paired, parametric):
    # get sample set
    sample_set = set(df[annot_colname])

    # insert each term into AnnotationHierarchy
    ah = AnnotationHierarchy(db, sample_set)

    for index, row in df.iterrows():
        id = row[annot_colname]
        if isinstance(db, NCBITaxonomyDb):
            id = int(id)
        intensity_list = row[samp_grps.all_intcols].tolist()
        ah.add_node(id, intensity_list)
    # define sample children and filter to informative
    ah.get_informative_nodes(min_peptides=min_peptides, min_children_non_leaf=min_children_non_leaf)
    intensity_all_ranks = ah.to_dataframe(samp_grps.all_intcols)

    results = common_stats(intensity_all_ranks, samp_grps, test, threshold, paired, parametric)
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