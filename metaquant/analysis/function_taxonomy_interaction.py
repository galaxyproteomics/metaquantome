from metaquant.databases.cog import cogCat
from metaquant.databases.cog import take_first_cog
from metaquant.util import stats, utils
import numpy as np
from metaquant.databases.NCBITaxonomyDb import NCBITaxonomyDb


def function_taxonomy_analysis(df, cog_name, lca_colname, samp_grps, threshold, data_dir):
    """
    Runs the function-taxonomy interaction analysis. The use of GO terms is required.
    The process is as follows:
    1. Reduce provided dataframe
    2. Normalize dataframe, so each GO term has its own row
    3. Map each GO term to its slim version - new column
    4. Drop duplicated peptide-GO slim combinations (so we don't double count)
    5. Drop rows with peptides that have an LCA taxon at a higher rank than the desired rank (`des_rank`)
    6. Map each taxon to its desired rank - new column
    7. Group by the new taxon column and the new GO term column
    :param df: joined taxonomy, intensity, and function tables
    :param cog_name: name of COG column in dataframe
    :param lca_colname: name of LCA column in dataframe.
    :param samp_grps: a SampleGroups object for this analysis
    :return: dataframe with taxon-function pairs and their associated total intensity
    """
    df = take_first_cog(df, cog_name)

    # select columns for adding purposes
    df_int = df[samp_grps.all_intcols + [cog_name] + [lca_colname]]

    # group by both cog and lca and add
    grouped = df_int.groupby(by=[cog_name, lca_colname]).sum(axis=1)

    # test
    # if test:
    #     results = stats.test_norm_intensity(grouped, samp_grps, threshold, paired, parametric)
    # else:
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
    if not data_dir:
        data_dir = utils.define_ontology_data_dir('taxonomy')
    ncbi = NCBITaxonomyDb(data_dir)
    results['rank'] = [ncbi.get_rank(int(elem)) for elem in taxids]
    results['taxon_name'] = ncbi.convert_taxid_to_name(taxids)

    return results

