
def split_func_list(golist, sep):
    """
    make set from a delimited string of GO terms
    :param golist: a string of GO terms separated by `sep`: ex. 'go1,go2,go3'
    :return: A set of GO terms as strings. Ex: {'go1', 'go2', 'go3'}
    """
    goset = set(golist.split(sep))
    return goset

def paste_func_set(goset, sep):
    """
    Write a set to a delimited string of GO terms
    :param goset: set of GOs
    :param sep: desired separator between GO terms in string
    :return: a string of GO terms, separated by `sep`
    """
    string = sep.join(goset)
    return string

def reduce_func(db, funclist, sep):
    split = split_func_list(funclist, sep)
    reduced = make_set_nonredundant(db, split)
    pasted = paste_func_set(reduced, sep)
    return pasted

def reduce_func_df(db, df, func_colname, sep):
    """
    Replace the old functional column with a nonredundant functional column
    :param df: Combined df
    :param func_colname: String, name of the functional column
    :param sep: String that separates the list of GO terms in the functional column
    :return: pandas DataFrame with nonredundant functional column
    """
    orig_func_col = df[func_colname]
    new_func_col = orig_func_col.apply(lambda x: reduce_func(db, x, sep))
    df[func_colname] = new_func_col
    return df

def make_set_nonredundant(db, funcset):
    """
    Take a set of go terms and return only the
    GO terms that are not the ancestor of any
    other term in the set
    :param goset: set of GO terms to be reduced
    :return: the set of GO terms with redundant terms removed
    """
    reduced_set = funcset.copy()
    for goid in funcset:
        ancestors = db.get_ancestors(goid)
        reduced_set.difference_update(ancestors)
    return reduced_set
