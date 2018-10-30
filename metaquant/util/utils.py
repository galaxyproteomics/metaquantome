import re
import pkg_resources
import os

DATA_DIR = pkg_resources.resource_filename('metaquant', 'data/')
GO_SUBDIR = 'go'
EC_SUBDIR = 'enzyme'
TAX_SUBDIR = 'ncbi'


def define_ontology_data_dir(ontology):
    base_ddir = DATA_DIR
    if ontology == "go":
        dir = GO_SUBDIR
    elif ontology == "ec":
        dir = EC_SUBDIR
    elif ontology == "taxonomy":
        dir = TAX_SUBDIR
    elif ontology == "cog":
        return base_ddir
    full_ddir = os.path.join(base_ddir, dir)
    return full_ddir


def safe_cast_to_list(obj):
    if not isinstance(obj, list):
        return [obj]
    else:
        return obj


def split_func_list(df, sep, func_colname):
    new_df = tidy_split(df, func_colname, sep=sep, keep=False)
    return new_df


def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row.
    modified from https://github.com/cognoma/genes/blob/721204091a96e55de6dcad165d6d8265e67e2a48/2.process.py#L61-L95

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df


def sniff_tax_names(df, tax_colname):  # todo: move to NCBI database
    '''
    if greater than 90% of entries contain numbers, then we say it is numeric
    otherwise, we attempt to convert the names to taxids
    :param df:
    :param tax_colname:
    :return:
    '''
    pattern = re.compile(r'[0-9]')  # little bit faster to compile
    is_numeric = df[tax_colname].str.contains(pattern)
    if is_numeric.mean() > 0.9:
        return False  # if any entries contain numbers, assume taxids (already converted missings to NA)
    else:
        return True  # else names
