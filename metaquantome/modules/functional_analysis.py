import metaquantome.modules.expand
import metaquantome.util.utils
from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
from metaquantome.databases.cog import cogCat
import metaquantome.databases.EnzymeDb as ec
from metaquantome.util import utils


def functional_analysis(df, func_colname, samp_grps, ontology, slim_down, data_dir):
    """
    Expand functional terms and aggregate intensities.

    :param df: A DataFrame. Missing values are 0. There may be multiple
    functional terms in each row.
    :param func_colname: The name of the column with functional terms.
    :param samp_grps: A SampleGroups() object.
    :param ontology: Functional ontology. Either 'go', 'ec', or 'cog'
    :param slim_down: Boolean. Whether to map terms to slim or not.
    :param data_dir: Directory to contain functional database files (ex. go-basic.obo)
    :return: A dataframe with a functional term and its associated sample-specific intensity in each
    row
    """
    db, norm_df = clean_function_df(data_dir, df, func_colname, ontology, slim_down)

    if ontology == "go":
        # filter to only those in full GO and non-missing
        df_clean = utils.filter_df(db, func_colname, norm_df)
        if slim_down:
            # map full GO terms to slims
            df_clean = slim_down_df(db, df_clean, func_colname)
        # expand
        results = metaquantome.modules.expand.\
            common_hierarchical_analysis(db, df_clean, func_colname, samp_grps)
        gos = [db.gofull[x] for x in results['id']]
        results['name'] = [x.name for x in gos]
        results['namespace'] = [x.namespace for x in gos]
    elif ontology == "ec":
        # filter to only those in Enzyme database and non-missing
        df_clean = utils.filter_df(db, func_colname, norm_df)
        results = metaquantome.modules.expand.\
            common_hierarchical_analysis(db, df_clean, func_colname, samp_grps)
        # use ec database to get term descriptions
        results['description'] = [db.ecdb[term]['descript'] for term in results.index]
    elif ontology == "cog":
        cog_sum_df = norm_df[[func_colname] + samp_grps.all_intcols].\
            groupby(func_colname).\
            sum()
        results = metaquantome.modules.expand.\
            common_hierarchical_analysis('cog', cog_sum_df, func_colname, samp_grps,
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
    # make a local copy, rather than modifying the original df
    df_loc = df_clean.copy()
    # get a panda Series object of the GO terms
    func_series = df_loc[go_colname]
    # convert the Series to a set, for map_set_to_slim
    func_set = set(func_series)
    # returns a dictionary with {<full_go_term>: <slim_go_term>, ... }
    mapper = godb.map_set_to_slim(func_set)
    # use map on the Series
    df_loc.loc[:, go_colname] = func_series.map(mapper)
    # return copy
    return df_loc


def clean_function_df(data_dir, df, func_colname, ontology, slim_down):
    """
    make functional terms nonredundant and normalize dataframe so there's only one functional
    term per row
    :param data_dir: directory to contain the database files for specified ontology
    :param df: DataFrame. May have multiple functional terms per row. Missing values
    should be 0
    :param func_colname: Name of column with functional terms
    :param ontology: Desired ontology.
    :param slim_down: Map full GO terms to the metagenomics slim. Applies for GO only.
    :return: A tuple of the database and the dataframe with one functional term in each row.
    """
    # db data dir
    if not data_dir:
        data_dir = utils.DATA_DIR
    # assign db
    db = None
    if ontology in {"go", "ec"}:
        if ontology == "go":
            db = GeneOntologyDb(data_dir, slim_down)
        elif ontology == "ec":
            db = ec.EnzymeDb(data_dir)
        # reduce df to non-redundant functional terms
        # todo: add sep to args (if desired at some point)
        red_df = metaquantome.util.utils.reduce_func_df(db=db, df=df, func_colname=func_colname, sep=',')
    else:
        red_df = df
    # normalize df, so each row has one functional term
    norm_df = utils.tidy_split(red_df, column=func_colname, sep=',')
    return db, norm_df
