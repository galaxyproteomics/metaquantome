import numpy as np

import metaquantome.analysis.expand
from metaquantome.util import utils
from metaquantome.analysis import functional_analysis as fa
from metaquantome.databases.NCBITaxonomyDb import NCBITaxonomyDb


def function_taxonomy_analysis(df, func_colname, pep_colname, ontology, overwrite, slim_down, tax_colname, samp_grps,
                               ft_tar_rank, ft_func_data_dir, ft_tax_data_dir):
    """
    Runs the function-taxonomy interaction analysis. The use of GO terms is required.
    The process is as follows:
    1. Reduce provided dataframe
    2. Normalize dataframe, so each GO term has its own row
    3. (optionally) Map each GO term to its slim version - new column
    4. Drop duplicated peptide-GO (slim) combinations (so we don't double count)
    5. Drop rows with peptides that have an LCA taxon at a higher rank than the desired rank (`des_rank`)
    6. Map each taxon to its desired rank - new column
    7. Group by the new taxon column and the new GO term column
    :param df: joined taxonomy, intensity, and function tables
    :param cog_name: name of COG column in dataframe
    :param tax_colname: name of LCA column in dataframe.
    :param samp_grps: a SampleGroups object for this analysis
    :return: dataframe with taxon-function pairs and their associated total intensity
    """
    # ---- arg checks ---- #
    # ontology must be go
    if ontology != 'go':
        ValueError('ontology must be "go" for function-taxonomy analysis')

    # ---- reduce, normalize, (optionally) slim ---- #
    godb, norm_df = fa.clean_function_df(ft_func_data_dir, df, func_colname, ontology, overwrite, slim_down)
    if slim_down:
        norm_df = fa.slim_down_df(godb, norm_df, func_colname)
    # remove peptide/go-term duplicates
    dedup_df = norm_df.\
        reset_index().\
        drop_duplicates(subset=[pep_colname, func_colname], keep='first').\
        set_index(pep_colname)
    # ---- get rank of lca ----- #
    # resolve data dir
    if not ft_tax_data_dir:
        ft_tax_data_dir = utils.define_ontology_data_dir('taxonomy')
    # load ncbi database
    ncbi = NCBITaxonomyDb(ft_tax_data_dir)
    # see if names. if so, convert to taxid
    if utils.sniff_tax_names(df, tax_colname):
        dedup_df[tax_colname] = ncbi.convert_name_to_taxid(dedup_df[tax_colname].tolist())
    else:
        dedup_df[tax_colname] = [int(x) for x in dedup_df[tax_colname]]
    # filter df to those that tax ids that non-NaN and are present in NCBI database
    is_not_nan = dedup_df[tax_colname].notnull()
    is_in_db = dedup_df[tax_colname].apply(ncbi.is_in_db)
    dedup_df = dedup_df.loc[is_not_nan & is_in_db]

    dedup_df['des_rank'] = dedup_df[tax_colname].apply(lambda x: des_rank_mapper(ft_tar_rank, x, ncbi))
    # filter out peptides that are less specific than query rank
    dedup_df = dedup_df[dedup_df['des_rank'] > 0]

    # ---- group by go and new des_rank column, then sum intensity ---- #
    # select columns for adding purposes
    df_int = dedup_df[samp_grps.all_intcols + [func_colname, 'des_rank']]
    # group by both cog and lca and add
    grouped = df_int.groupby(by=[func_colname, 'des_rank']).sum(axis=1)
    # get groupwise counts (i.e., unique peptides)
    counts = df_int.groupby(by=[func_colname, 'des_rank']).size().to_frame(name="n_peptide")
    ints_and_counts = grouped.join(counts)

    # ---- output prep ---- #
    # calculate group means
    results = metaquantome.analysis.expand.calc_means(ints_and_counts, samp_grps)
    # take log of intensities for return
    results[results == 0] = np.nan
    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])
    # split the cog/lca index back into 2 columns
    results.reset_index(inplace=True)
    # add go description
    gos = [godb.gofull[x] for x in results[func_colname]]
    results['name'] = [x.name for x in gos]
    results['namespace'] = [x.namespace for x in gos]
    # translate ids back to names
    taxids = results['des_rank']
    # get ranks
    results['rank'] = [ncbi.get_rank(int(elem)) for elem in taxids]
    results['taxon_name'] = ncbi.convert_taxid_to_name(taxids)
    return results


def des_rank_mapper(des_rank, taxid, ncbi):
    """
    function for mapping a taxid to a desired rank.
    if des_rank is lower than rank of taxid, 0 is returned (to be filtered out later)
    :param des_rank:
    :param taxid:
    :param ncbi:
    :return:
    """
    dict_mapper = ncbi.map_id_to_desired_ranks([des_rank], int(taxid))
    if len(dict_mapper) == 1:
        return dict_mapper[des_rank]
    else:
        return 0


