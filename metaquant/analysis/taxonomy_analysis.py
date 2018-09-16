from metaquant.databases.NCBITaxonomyDb import NCBITaxonomyDb
import metaquant.analysis.common as cha
from metaquant.util import utils


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
    # # get sample set
    # sample_set = set(df[tax_colname])
    #
    # # insert each term into AnnotationHierarchy
    # ah = AnnotationHierarchy(ncbi, sample_set)
    # for index, row in df.iterrows():
    #     taxid_int = int(row[tax_colname])
    #     intensity_list = row[samp_grps.all_intcols].tolist()
    #     ah.add_node(taxid_int, intensity_list)
    #
    # # define sample children and filter to informative
    # ah.get_informative_nodes(min_peptides=min_peptides, min_children_non_leaf=min_children_non_leaf)
    #
    # intensity_all_ranks = ah.to_dataframe(samp_grps.all_intcols)
    #
    # # test
    # if test and samp_grps.ngrps == 2:
    #     results = stats.test_norm_intensity(intensity_all_ranks, samp_grps, threshold, paired, parametric)
    # else:
    #     results = stats.calc_means(intensity_all_ranks, samp_grps)
    #
    # # replace nan with zero, so that np.log2 returns nan
    # results[results == 0] = np.nan
    #
    # # take log of intensities for return
    # results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])
    results['rank'] = results['id'].apply(ncbi.get_rank)

    # translate ids back to names
    results['taxon_name'] = ncbi.convert_taxid_to_name(results['id'])

    return results
