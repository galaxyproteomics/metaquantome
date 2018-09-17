from metaquant.databases.NCBITaxonomyDb import NCBITaxonomyDb
import metaquant.analysis.common as cha
from metaquant.util import utils
import pandas as pd


def taxonomy_analysis(df, samp_grps, test, threshold, paired, parametric, data_dir, tax_colname='lca',
                      min_peptides=0,
                      min_children_non_leaf=0):
    if not data_dir:
        data_dir = utils.define_ontology_data_dir('taxonomy')
    # load ncbi database
    ncbi = NCBITaxonomyDb(data_dir)

    # check for numeric characters, which indicates taxid
    # if is name, convert to taxid
    # keep as character until querying ncbi database
    if utils.sniff_tax_names(df, tax_colname):
        df[tax_colname] = ncbi.convert_name_to_taxid(df[tax_colname])

    results = cha.common_hierarchical_analysis(ncbi, df, tax_colname, samp_grps,
                                               min_peptides, min_children_non_leaf,
                                               test, threshold, paired, parametric)
    results['rank'] = results['id'].apply(ncbi.get_rank)

    # translate ids back to names
    results['taxon_name'] = ncbi.convert_taxid_to_name(results['id'])
    return results
