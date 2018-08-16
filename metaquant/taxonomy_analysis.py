import pandas as pd
from metaquant import stats
from metaquant import taxonomy_database
import numpy as np

from metaquant.taxonomy_database import BASIC_TAXONOMY_TREE


def taxonomy_analysis(df, samp_grps, test, threshold, paired, parametric, data_dir, tax_colname='lca', overwrite=False):

    # load ncbi database
    ncbi = taxonomy_database.ncbi_database_handler(data_dir)

    # taxids, uniqify
    lca = df[tax_colname].unique()

    # get full lineage for each unique taxid
    lineages = [pd.DataFrame(
        taxonomy_database.map_id_to_desired_ranks(BASIC_TAXONOMY_TREE, taxid, ncbi),
        index=[taxid]) for taxid in lca
    ]
    full_lineage = pd.concat(lineages, sort=False)

    # map lineages to df with intensity
    joined = df.join(full_lineage, on=[tax_colname])

    # new - just add up through ranks
    intensity_all_ranks = pd.concat([stats.group_and_sum_by_rank(joined, x, samp_grps.all_intcols, norm_to_rank=False) for x in
                                     BASIC_TAXONOMY_TREE])

    # test
    if test and samp_grps.ngrps == 2:
        results = stats.test_norm_intensity(intensity_all_ranks, samp_grps, threshold, paired, parametric)
    else:
        results = stats.calc_means(intensity_all_ranks, samp_grps)

    # replace nan with zero, so that np.log2 returns nan
    results[results == 0] = np.nan

    # take log of intensities for return
    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])

    # translate ids back to names
    results['taxon_name'] = taxonomy_database.convert_taxid_to_name(results['id'], ncbi)

    return results
