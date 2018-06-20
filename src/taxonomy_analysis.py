import pandas as pd
from src import stats
from src import phylo_tree
import numpy as np


BASIC_TAXONOMY_TREE = ["phylum",
                       "class",
                       "order",
                       "family",
                       "genus",
                       "species"]


def taxonomy_analysis(df, samp_grps, test, threshold, paired, tax_colname='lca'):

    # load ncbi database
    ncbi = phylo_tree.load_ncbi()

    # taxids, uniqify
    lca = df[tax_colname].unique()

    # get full lineage for each unique taxid
    lineages = [pd.DataFrame(
        phylo_tree.get_desired_ranks_from_lineage(BASIC_TAXONOMY_TREE, taxid, ncbi),
        index=[taxid]) for taxid in lca
    ]
    full_lineage = pd.concat(lineages)

    # map lineages to df with intensity
    joined = df.join(full_lineage, on=[tax_colname])

    # new - just add up through ranks
    intensity_all_ranks = pd.concat([abundance_rank(joined, x, samp_grps.all_intcols, norm_to_rank=False) for x in BASIC_TAXONOMY_TREE])

    # test
    if test and samp_grps.ngrps == 2:
        results = stats.test_norm_intensity(intensity_all_ranks, samp_grps, threshold, paired)
    else:
        results = stats.calc_means(intensity_all_ranks, samp_grps)

    # replace nan with zero, so that np.log2 returns nan
    results[results == 0] = np.nan

    # take log of intensities for return
    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])

    # translate ids back to names
    results['taxon_name'] = phylo_tree.convert_taxid_to_name(results['id'], ncbi)

    return results


def abundance_rank(df, rank, all_intcols, norm_to_rank=False):
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
