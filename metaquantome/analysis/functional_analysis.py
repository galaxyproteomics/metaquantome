import metaquantome.analysis.expand
from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
from metaquantome.databases.cog import cogCat
from metaquantome.databases.cog import take_first_cog
import metaquantome.databases.EnzymeDb as ec
from metaquantome.util import utils, funcutils


def functional_analysis(df, func_colname, samp_grps, ontology, slim_down, data_dir, overwrite):
    """
    Expand functional terms and aggregate intensities.
    :param df: A DataFrame. Missing values are 0. There may be multiple
    functional terms in each row.
    :param func_colname: The name of the column with functional terms.
    :param samp_grps: A SampleGroups() object.
    :param ontology: Functional ontology. Either 'go', 'ec', or 'cog'
    :param slim_down: Boolean. Whether to map terms to slim or not.
    :param data_dir: Directory to contain functional database files (ex. go-basic.obo)
    :param overwrite: Whether to update the database files (applies to EC and GO only)
    :return: A dataframe with a functional term and its associated sample-specific intensity in each
    row
    """
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


def slim_down_df(godb, df_clean, go_colname):
    """
    Maps GO column to slim GO
    :param godb: The GO database
    :param df_clean: DataFrame with one GO term per row
    :param go_colname: Name for the column with GO terms
    :return: DataFrame with old GO terms replaced with slim GO terms
    """
    df_loc = df_clean.copy()
    func_series = df_loc[go_colname]
    func_set = set(func_series)
    mapper = godb.map_set_to_slim(func_set)
    df_loc.loc[:, go_colname] = func_series.map(mapper)
    return df_loc


def filter_func_df(db, func_colname, norm_df):
    """
    Filter DataFrame to non-missing and valid terms
    :param db: Relevant database
    :param func_colname: Name of column with functional terms
    :param norm_df: Dataframe with one functional term per row
    :return: DataFrame with only non-missing terms and those present in the database
    """
    is_not_nan = ~norm_df[func_colname].isnull()
    is_in_db = norm_df[func_colname].apply(db.is_in_db)
    df_clean = norm_df.loc[is_not_nan & is_in_db].copy(deep=True)  # copy() avoids setting with copy warning
    return df_clean


def clean_function_df(data_dir, df, func_colname, ontology, overwrite, slim_down):
    """
    make functional terms nonredundant and normalize dataframe so there's only one functional
    term per row
    :param data_dir: directory to contain the database files for specified ontology
    :param df: DataFrame. May have multiple functional terms per row. Missing values
    should be 0
    :param func_colname: Name of column with functional terms
    :param ontology: Desired ontology.
    :param overwrite: Whether to update the database.
    :param slim_down: Map full GO terms to the metagenomics slim. Applies for GO only.
    :return: A tuple of the database and the dataframe with one fuctional term in each row.
    """
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


