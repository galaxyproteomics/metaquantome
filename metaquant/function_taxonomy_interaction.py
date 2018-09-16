from metaquant.cog import cogCat
from metaquant.cog import take_first_cog
from metaquant import stats
import numpy as np
from metaquant.NCBITaxonomyDb import NCBITaxonomyDb
from metaquant import utils


def function_taxonomy_analysis(df, cog_name, lca_colname, samp_grps, test, threshold, paired, parametric, data_dir):
    """
    Runs the function-taxonomy interaction analysis. For the documentation of other arguments, see metaquant.py
    :param df: joined taxonomy, intensity, and function tables
    :param cog_name: name of COG column in dataframe
    :param lca_colname: name of LCA column in dataframe.
    :param samp_grps: a SampleGroups object for this analysis
    :return: dataframe with taxon-function pairs and their associated total intensity
    """
    # todo: add option for lca or rank-level
    # take first cog
    df = take_first_cog(df, cog_name)

    # select columns for adding purposes
    df_int = df[samp_grps.all_intcols + [cog_name] + [lca_colname]]

    # group by both cog and lca and add
    grouped = df_int.groupby(by=[cog_name, lca_colname]).sum(axis=1)

    # test
    if test:
        results = stats.test_norm_intensity(grouped, samp_grps, threshold, paired, parametric)
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
    if not data_dir:
        data_dir = utils.define_ontology_data_dir('taxonomy')
    ncbi = NCBITaxonomyDb(data_dir)
    results['rank'] = [ncbi.get_rank(int(elem)) for elem in taxids]
    results['taxon_name'] = ncbi.convert_taxid_to_name(taxids)

    return results

