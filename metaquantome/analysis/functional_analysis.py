import metaquantome.analysis.expand
from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
from metaquantome.databases.cog import cogCat
from metaquantome.databases.cog import take_first_cog
import metaquantome.databases.EnzymeDb as ec
from metaquantome.util import utils, funcutils


def functional_analysis(df, func_colname, samp_grps, ontology, slim_down, data_dir, overwrite):
    db, norm_df = clean_function_df(data_dir, df, func_colname, ontology, overwrite, slim_down)

    if ontology == "go":
        # filter to only those in full GO and non-missing
        df_clean = filter_func_df(db, func_colname, norm_df)
        if slim_down:
            df_clean = slim_down_df(db, df_clean, func_colname)
        results = metaquantome.analysis.expand.common_hierarchical_analysis(db, df_clean, func_colname, samp_grps)
        # todo: replace hard coded column names in utils
        gos = [db.gofull[x] for x in results['id']]
        results['name'] = [x.name for x in gos]
        results['namespace'] = [x.namespace for x in gos]
    elif ontology == "ec":
        # filter to only those in Enzyme database and non-missing
        df_clean = filter_func_df(db, func_colname, norm_df)
        results = metaquantome.analysis.expand.common_hierarchical_analysis(db, df_clean, func_colname, samp_grps)
        results['description'] = [db.ecdb[term]['descript'] for term in results.index]
    elif ontology == "cog":
        cog_df = take_first_cog(df, func_colname)
        cog_sum_df = cog_df[[func_colname] + samp_grps.all_intcols].\
            groupby(func_colname).\
            sum()
        results = metaquantome.analysis.expand.common_hierarchical_analysis('cog', cog_sum_df, func_colname, samp_grps,
                                                                            hierarchical=False)
        results['description'] = [cogCat[x] for x in results.index]
    else:
        raise ValueError("the desired ontology is not supported. " +
                         "Please use either GO (ontology = 'go'), " +
                         "COG (ontology = 'cog'), or EC numbers (ontology = 'ec')")
    return results


def slim_down_df(db, df_clean, func_colname):
    df_loc = df_clean.copy()
    func_series = df_loc[func_colname]
    func_set = set(func_series)
    mapper = db.map_set_to_slim(func_set)
    df_loc.loc[:, func_colname] = func_series.map(mapper)
    return df_loc


def filter_func_df(db, func_colname, norm_df):
    is_not_nan = ~norm_df[func_colname].isnull()
    is_in_db = norm_df[func_colname].apply(db.is_in_db)
    df_clean = norm_df.loc[is_not_nan & is_in_db].copy(deep=True)  # copy() avoids setting with copy warning
    return df_clean


def clean_function_df(data_dir, df, func_colname, ontology, overwrite, slim_down):
    # db data dir
    if not data_dir:
        data_dir = utils.define_ontology_data_dir(ontology)
    # assign db
    db = None
    if ontology in {"go", "ec"}:
        if ontology == "go":
            db = GeneOntologyDb(data_dir, slim_down, overwrite)
        elif ontology == "ec":
            db = ec.EnzymeDb(data_dir, overwrite)
        # reduce df to non-redundant functional terms
        # todo: add sep to args
        red_df = funcutils.reduce_func_df(db=db, df=df, func_colname=func_colname, sep=',')
    else:
        red_df = df
    # normalize df, so each row has one functional term
    norm_df = utils.split_func_list(red_df, sep=',', func_colname=func_colname)
    return db, norm_df


