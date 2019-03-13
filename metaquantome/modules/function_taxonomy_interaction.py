import numpy as np

import metaquantome.modules.expand
from metaquantome.util import utils
from metaquantome.modules import functional_analysis as fa
from metaquantome.databases.NCBITaxonomyDb import NCBITaxonomyDb


def function_taxonomy_analysis(df, func_colname, pep_colname, ontology, slim_down, tax_colname, samp_grps, ft_tar_rank,
                               data_dir):
    # todo: add normalization module for ft modules
    # choose bp, cc, or mf
    # don't return bp, cc, or mf themselves
    """
    Runs the function-taxonomy interaction modules. The use of GO terms is required.
    The process is as follows:
    1. Reduce provided dataframe
    2. Normalize dataframe, so each GO term has its own row
    3. (optionally) Map each GO term to its slim version - new column
    4. Drop duplicated peptide-GO (slim) combinations (so we don't double count)
    5. Drop rows with peptides that have an LCA taxon at a higher rank than the desired rank (`des_rank`)
    6. Map each taxon to its desired rank - new column
    7. Group by the new taxon column and the new GO term column
    :param df: joined taxonomy, intensity, and function tables
    :param func_colname: name of function column in dataframe
    :param pep_colname: name of peptide column in dataframe
    :param ontology: name of functional ontology. must be 'go'
    :param slim_down: whether to map full GO terms to metagenomics slim GO terms
    :param tax_colname: name of LCA column in dataframe.
    :param samp_grps: a SampleGroups object for this modules
    :param ft_tar_rank: rank at which to group taxonomy. Default is 'genus'
    :param data_dir: data directory
    :return: dataframe with taxon-function pairs and their associated total intensity
    """
    # ---- arg checks ---- #
    # ontology must be go
    if ontology != 'go':
        ValueError('ontology must be "go" for function-taxonomy modules')

    # ---- reduce, normalize, (optionally) slim ---- #
    godb, norm_df = fa.clean_function_df(data_dir, df, func_colname, ontology, slim_down)
    if slim_down:
        norm_df = fa.slim_down_df(godb, norm_df, func_colname)
    # remove peptide/go-term duplicates (in the case that different GO term annotations
    # for the same peptide are mapped to the same slim GO term)
    # index is named peptide
    dedup_df = norm_df.\
        reset_index().\
        drop_duplicates(subset=['peptide', func_colname], keep='first').\
        set_index('peptide')
    # ---- get rank of lca ----- #
    # resolve data dir
    if not data_dir:
        data_dir = utils.DATA_DIR
    # load ncbi database
    ncbi = NCBITaxonomyDb(data_dir)
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
    # filter out peptides that are less specific than query rank (which have a taxid of 0)
    dedup_df = dedup_df[dedup_df['des_rank'] > 0]

    # ---- group by go and new des_rank column, then sum intensity ---- #
    # select columns for adding purposes
    df_int = dedup_df[samp_grps.all_intcols + [func_colname, 'des_rank']]
    # do counts
    df_counts = df_int.copy()
    df_counts.loc[:, samp_grps.all_intcols] = df_counts.loc[:, samp_grps.all_intcols] > 0
    # group by both cog and lca and add
    grouped = df_int.groupby(by=[func_colname, 'des_rank']).sum(axis=1)
    # get groupwise counts (i.e., unique peptides)
    # multiply by 1 to convert any single booleans (True) to 1
    counts = df_counts.groupby(by=[func_colname, 'des_rank']).sum(axis=1) * 1
    ints_and_counts = grouped.join(counts, rsuffix='_n_peptide')

    # ---- output prep ---- #
    # calculate group means
    results = metaquantome.modules.expand.calc_means(ints_and_counts, samp_grps)
    # take log of intensities for return
    results[results == 0] = np.nan
    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])
    # split the cog/lca index back into 2 columns
    results.reset_index(inplace=True)
    # add go description
    gos = [godb.gofull[x] for x in results[func_colname]]
    results['go_id'] = results[func_colname]
    results['name'] = [x.name for x in gos]
    results['namespace'] = [x.namespace for x in gos]
    # translate ids back to names
    taxids = results['des_rank']
    # get ranks
    results['tax_id'] = taxids
    results['rank'] = [ncbi.get_rank(int(elem)) for elem in taxids]
    results['taxon_name'] = ncbi.convert_taxid_to_name(taxids)
    # drop des_rank column
    results.drop('des_rank', axis=1, inplace=True)
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
    dict_mapper = ncbi.map_id_to_desired_ranks([des_rank], int(taxid))  # todo: fix if taxid is already int
    if len(dict_mapper) == 1:  # this should always just be 1 or 0, because we're querying with one term
        return dict_mapper[des_rank]
    else:
        return 0
