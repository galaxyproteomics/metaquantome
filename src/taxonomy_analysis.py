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


def taxonomy_analysis(df, samp_grps, test, threshold, paired):

    # load ncbi database
    ncbi = phylo_tree.load_ncbi()

    # get lineage of lca, make df
    lca = df['lca']

    full_lineage = pd.concat(
        [pd.DataFrame(
            phylo_tree.get_desired_ranks_from_lineage(BASIC_TAXONOMY_TREE, taxid, ncbi),
            index=[taxid]) for taxid in df['lca']
        ]
    )

    # join to df with intensities
    joined = df.join(full_lineage, on='lca')

    # old - norm by rank
    # norm_intensity_all_ranks = pd.concat([rel_abundance_rank(joined, x, samp_grps.all_intcols) for x in BASIC_TAXONOMY_TREE])

    # new - just add up through ranks
    intensity_all_ranks = pd.concat([abundance_rank(joined, x, samp_grps.all_intcols) for x in BASIC_TAXONOMY_TREE])

    # test
    if test:
        results = stats.test_norm_intensity(intensity_all_ranks, samp_grps, threshold, paired)
    else:
        results = stats.calc_means(intensity_all_ranks, samp_grps)

    # take log of intensities for return
    # replace zeros with nan
    results[results == 0] = np.nan

    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])

    # translate ids back to names
    results['id'] = phylo_tree.convert_taxid_to_name(results['id'], ncbi)

    return results


def abundance_rank(df, rank, all_intcols, norm_to_rank=False):
    summed_abund = df.groupby(by=rank)[all_intcols].sum(axis=0)

    # normalize to each sample
    if norm_to_rank:
        return_df = summed_abund / summed_abund.sum(axis=0)
    else:
        return_df = summed_abund

    return_df['rank'] = rank
    return_df['id'] = return_df.index
    return return_df
