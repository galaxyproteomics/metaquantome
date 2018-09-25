from metaquant.databases.GeneOntologyDb import GeneOntologyDb
from metaquant.databases.cog import cogCat
from metaquant.databases.cog import take_first_cog
import metaquant.databases.EnzymeDb as ec
from metaquant.util import utils
import metaquant.analysis.common as cha


def functional_analysis(df, func_colname, samp_grps, test, threshold, ontology, slim_down, paired, parametric, data_dir,
                        overwrite, min_peptides=0, min_children_non_leaf=0):

    norm_df = utils.split_func_list(df, sep=',', func_colname=func_colname)
    if not data_dir:
        data_dir = utils.define_ontology_data_dir(ontology)
    if ontology == "go":
        go = GeneOntologyDb(data_dir, slim_down, overwrite)
        if slim_down:
            func_series = norm_df[func_colname]
            func_set = set(func_series)
            mapper = go.map_set_to_slim(func_set)
            norm_df[func_colname] = norm_df[func_colname].map(mapper)
        results = cha.common_hierarchical_analysis(go, norm_df, func_colname, samp_grps,
                                                   min_peptides, min_children_non_leaf,
                                                   test, threshold, paired, parametric)
        # todo: replace hard coded column names in utils
        gos = [go._safe_query_go(x) for x in results['id']]
        results['name'] = [x.name for x in gos]
        results['namespace'] = [x.namespace for x in gos]
    elif ontology == "cog":
        cog_df = take_first_cog(df, func_colname)
        cog_sum_df = cog_df[[func_colname] + samp_grps.all_intcols].\
            groupby(func_colname).\
            sum()
        df_to_return = cog_sum_df
        df_to_return['description'] = [cogCat[x] for x in df_to_return.index]
        results = cha.common_stats(df_to_return, samp_grps, test, threshold, paired, parametric)
    elif ontology == "ec":
        ec_db = ec.EnzymeDb(data_dir, overwrite)
        results = cha.common_hierarchical_analysis(ec_db, norm_df, func_colname, samp_grps,
                                                   min_peptides, min_children_non_leaf,
                                                   test, threshold, paired, parametric)
        results['descript'] = [ec_db.ecdb[term]['descript'] for term in results.index]
    else:
        raise ValueError("the desired ontology is not supported. " +
                         "Please use either GO (ontology = 'go'), " +
                         "COG (ontology = 'cog'), or EC numbers (ontology = 'ec')")
    return results


