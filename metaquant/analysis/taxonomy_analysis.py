from metaquant.databases.NCBITaxonomyDb import NCBITaxonomyDb
import metaquant.analysis.common as cha
from metaquant.util import utils


def taxonomy_analysis(df, samp_grps, data_dir, tax_colname='lca', min_peptides=0, min_children_non_leaf=0, threshold=0):
    if not data_dir:
        data_dir = utils.define_ontology_data_dir('taxonomy')
    # load ncbi database
    ncbi = NCBITaxonomyDb(data_dir)

    # check for numeric characters, which indicates taxid
    # if is name, convert to taxid
    # keep as character until querying ncbi database
    if utils.sniff_tax_names(df, tax_colname):
        df[tax_colname] = ncbi.convert_name_to_taxid(df[tax_colname])
    else:
        df[tax_colname] = [int(x) for x in df[tax_colname]]
    # filter df to those that tax ids that non-NaN and are present in NCBI database
    is_not_nan = df[tax_colname].notnull()
    is_in_db = df[tax_colname].apply(ncbi.is_in_db)
    df_clean = df.loc[is_not_nan & is_in_db]
    results = cha.common_hierarchical_analysis(ncbi, df_clean, tax_colname, samp_grps, min_peptides,
                                               min_children_non_leaf, threshold)
    results['rank'] = results['id'].apply(ncbi.get_rank)

    # translate ids back to names
    results['taxon_name'] = ncbi.convert_taxid_to_name(results['id'])
    return results
