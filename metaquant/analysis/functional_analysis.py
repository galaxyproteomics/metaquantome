from metaquant.databases.GeneOntologyDb import GeneOntologyDb
from metaquant.databases.cog import cogCat
from metaquant.databases.cog import take_first_cog
import metaquant.databases.EnzymeDb as ec
from metaquant.util import utils
import metaquant.analysis.common as cha


def functional_analysis(df, func_colname, samp_grps, ontology, slim_down, data_dir, overwrite,
                        min_peptides=0, min_children_non_leaf=0, threshold=0):
    norm_df = utils.split_func_list(df, sep=',', func_colname=func_colname)
    if not data_dir:
        data_dir = utils.define_ontology_data_dir(ontology)
    if ontology == "go":
        go = GeneOntologyDb(data_dir, slim_down, overwrite)
        # filter to only those in full GO
        is_not_nan = ~norm_df[func_colname].isnull()
        is_in_db = norm_df[func_colname].apply(go.is_in_db)
        df_clean = norm_df.loc[is_not_nan & is_in_db].copy(deep=True)  # copy() avoids setting with copy warning
        if slim_down:
            func_series = df_clean[func_colname]
            func_set = set(func_series)
            mapper = go.map_set_to_slim(func_set)
            df_clean.loc[:, func_colname] = func_series.map(mapper)
        results = cha.common_hierarchical_analysis(go, df_clean, func_colname, samp_grps,
                                                   min_peptides, min_children_non_leaf,
                                                   threshold)
        # todo: replace hard coded column names in utils
        gos = [go.gofull[x] for x in results['id']]
        results['name'] = [x.name for x in gos]
        results['namespace'] = [x.namespace for x in gos]

    elif ontology == "cog":
        # todo: change to common hierarchical analysis
        cog_df = take_first_cog(df, func_colname)
        cog_sum_df = cog_df[[func_colname] + samp_grps.all_intcols].\
            groupby(func_colname).\
            sum()
        df_to_return = cog_sum_df
        df_to_return['description'] = [cogCat[x] for x in df_to_return.index]
        results = df_to_return

    elif ontology == "ec":
        ec_db = ec.EnzymeDb(data_dir, overwrite)
        # filter to only those in Enzyme database
        is_not_nan = ~norm_df[func_colname].isnull()
        is_in_db = norm_df[func_colname].apply(ec_db.is_in_db)
        df_clean = norm_df.loc[is_not_nan & is_in_db].copy(deep=True)  # copy() avoids setting with copy warning
        results = cha.common_hierarchical_analysis(ec_db, df_clean, func_colname, samp_grps,
                                                   min_peptides, min_children_non_leaf,
                                                   threshold)
        results['descript'] = [ec_db.ecdb[term]['descript'] for term in results.index]
    else:
        raise ValueError("the desired ontology is not supported. " +
                         "Please use either GO (ontology = 'go'), " +
                         "COG (ontology = 'cog'), or EC numbers (ontology = 'ec')")
    return results


