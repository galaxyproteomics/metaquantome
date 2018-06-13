import pandas as pd
from src import stats
from src import phylo_tree

FULL_TAXONOMIC_TREE = ["superkingdom",
                       "kingdom",
                       "subkingdom",
                       "superphylum",
                       "phylum",
                       "subphylum",
                       "superclass",
                       "class",
                       "subclass"
                       "infraclass",
                       "superorder",
                       "order",
                       "suborder",
                       "infraorder",
                       "parvorder",
                       "superfamily",
                       "family",
                       "subfamily",
                       "tribe",
                       "subtribe",
                       "genus",
                       "subgenus",
                       "species_group",
                       "species_subgroup",
                       "species",
                       "subspecies",
                       "varietas",
                       "forma"]

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
        [pd.DataFrame(phylo_tree.get_desired_ranks_from_lineage(BASIC_TAXONOMY_TREE, taxid, ncbi), index=[taxid]) for taxid in df['lca']]
    )

    # join to df with intensities
    joined = df.join(full_lineage, on='lca')

    # add up through ranks
    norm_intensity_all_ranks = pd.concat([rel_abundance_rank(joined, x, samp_grps.all_intcols) for x in BASIC_TAXONOMY_TREE])

    # test
    if test:
        results = stats.test_norm_intensity(norm_intensity_all_ranks, samp_grps, threshold, paired)
    else:
        results = stats.calc_means(norm_intensity_all_ranks, samp_grps)

    # translate ids back to names
    results['id'] = phylo_tree.convert_taxid_to_name(results['id'], ncbi)

    return results


def rel_abundance_rank(df, rank, all_intcols):
    summed_abund = df.groupby(by=rank)[all_intcols].sum(axis=0)

    # normalize to each sample
    rel_abundance = summed_abund / summed_abund.sum(axis=0)
    rel_abundance['rank'] = rank
    rel_abundance['id'] = rel_abundance.index
    return rel_abundance
