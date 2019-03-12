import re
import pkg_resources
import os
from urllib import request

BASE_DIR = pkg_resources.resource_filename('metaquantome', '/')
DATA_DIR = os.path.join(BASE_DIR, 'data')
MISSING_VALUES = ["", "0", "NA", "NaN", "0.0"]
ONTOLOGIES = ['cog', 'go', 'ec']
P_COLNAME = 'p'
P_CORR_COLNAME = 'corrected_p'
TEST_DIR = os.path.join(DATA_DIR, 'test')


def stream_to_file_from_url(url, tar):
    """
    uses the urllib library to open a URL and write the contents
    :param url: url of remote file to read
    :param tar: target local file
    :return: None
    """
    f = request.urlopen(url)
    data = f.read()
    with open(tar, 'wb') as tarfile:
        tarfile.write(data)
    f.close()


def safe_cast_to_list(string):
    """
    this is used to cast a string object to a length-one list,
    which is needed for some pandas operations.
    it also detects if something is already a list,
    and then does nothing (as to avoid a double list)

    :param string: string
    :return: either the string itself or a list comprising the string
    """
    if not isinstance(string, list):
        return [string]
    else:
        return string


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
    """
    if greater than 90% of entries contain numbers, then we say it is numeric
    otherwise, we attempt to convert the names to taxids
    :param df:
    :param tax_colname:
    :return:
    """
    pattern = re.compile(r'[0-9]')  # little bit faster to compile
    is_numeric = df[tax_colname].str.contains(pattern)
    if is_numeric.mean() > 0.9:
        return False  # if any entries contain numbers, assume taxids (already converted missings to NA)
    else:
        return True  # else names


def filter_df(db, annot_colname, norm_df):
    """
    Filter DataFrame to non-missing and valid terms

    :param db: Relevant database
    :param annot_colname: Name of column with terms (NCBI, GO, or EC)
    :param norm_df: Dataframe with one term per row
    :return: DataFrame with only non-missing terms and those present in the database
    """
    is_not_nan = ~norm_df[annot_colname].isnull()
    is_in_db = norm_df[annot_colname].apply(db.is_in_db)
    df_clean = norm_df.loc[is_not_nan & is_in_db].copy(deep=True)  # copy() avoids setting with copy warning
    return df_clean


def reduce_func(db, funclist, sep):
    """
    Take a set of functional terms and return only the
     terms that are not the ancestor of any
    other term in the set.

    :param db: reference database
    :param funclist: string that contains functional terms separated by <sep>
    :param sep: character separating the terms in the list
    :return:
    """
    split = set(funclist.split(sep))
    reduced = split.copy()
    for goid in split:
        ancestors = db.get_ancestors(goid)
        reduced.difference_update(ancestors)
    pasted = sep.join(reduced)
    return pasted


def reduce_func_df(db, df, func_colname, sep):
    """
    Replace the old functional column with a nonredundant functional column

    :param df: Combined df
    :param func_colname: String, name of the functional column
    :param sep: String that separates the list of terms in the functional column
    :return: pandas DataFrame with nonredundant functional column
    """
    orig_func_col = df[func_colname]
    new_func_col = orig_func_col.apply(lambda x: reduce_func(db, x, sep))
    df[func_colname] = new_func_col
    return df
