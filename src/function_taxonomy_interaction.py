from src.cog import cogCat
from src.cog import take_first_cog
from src import stats
import numpy as np
from src import phylo_tree


def function_taxonomy_analysis(df, cog_name, lca_colname, samp_grps, test, threshold, paired):

    # take first cog
    df = take_first_cog(df, cog_name)

    # select columns for adding purposes
    df_int = df[samp_grps.all_intcols + [cog_name] + [lca_colname]]

    # group by both cog and lca and add
    grouped = df_int.groupby(by=[cog_name, lca_colname]).sum(axis=1)

    # test
    if test:
        results = stats.test_norm_intensity(grouped, samp_grps, threshold, paired)
    else:
        results = stats.calc_means(grouped, samp_grps)

    # take log of intensities for return
    results[results == 0] = np.nan
    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])

    # split the cog/lca index back into 2 columns
    results.reset_index(inplace=True)

    # add cog description
    results['cog_descript'] = [cogCat[x] for x in results[cog_name]]

    # translate ids back to names
    taxids = results[lca_colname]

    # get ranks
    ncbi = phylo_tree.load_ncbi()
    taxid_to_rank_dict = ncbi.get_rank(taxids)
    results['rank'] = [taxid_to_rank_dict[int(elem)] for elem in taxids]

    results['taxon_name'] = phylo_tree.convert_taxid_to_name(taxids, ncbi)

    return results

