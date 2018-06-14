import pandas as pd
import itertools
from src import phylo_tree
import re
import numpy as np
from io import StringIO


MISSING_VALUES = ["", "0", "NA", "NaN", "0.0"]
ONTOLOGIES = ['cog', 'go']


def read_intensity_table(file, samp_grps, pep_colname):
    # read in data
    df = pd.read_table(file, sep="\t", index_col=pep_colname,
                       dtype=samp_grps.dict_numeric_cols,
                       na_values=MISSING_VALUES,
                       low_memory=False)

    # drop columns where all are NA
    df.dropna(axis=1, how="all", inplace=True)

    # change missing intensities to 0
    values = {x: 0 for x in samp_grps.all_intcols}
    df.fillna(values, inplace=True)

    # then change missing taxa/function to 'unknown'
    df.fillna('unknown', inplace=True)

    return df


def read_taxonomy_table(file, pep_colname, tax_colname):
    # always read as character
    df = pd.read_table(file, sep="\t", index_col=pep_colname,
                       na_values=MISSING_VALUES, dtype={'lca': object})

    # check for numeric characters, which indicates taxid
    # if numeric, convert
    if sniff_tax_names(df):
        ncbi = phylo_tree.load_ncbi()
        df['lca'] = phylo_tree.convert_name_to_taxid(df['lca'], ncbi)

    return df


def sniff_tax_names(df):
    pattern = re.compile(r'[0-9]')  # little bit faster to compile
    is_numeric = df['lca'].str.contains(pattern)
    if is_numeric.any():
        return False # if any entries contain numbers, assume taxids (already converted missings to NA)
    else:
        return True # else names


def read_function_table(file, pep_colname, ontology):
    """

    :param file:
    :param pep_colname:
    :param func_names: must be *list* of function names
    :return:
    """
    df = pd.read_table(file, sep="\t", index_col=pep_colname,
                       na_values=MISSING_VALUES)

    # make sure that ontology is supported
    if ontology not in set(df):
        raise ValueError('function columns do not have correct names, "go" and/or "cog"')

    return(df)


def join_on_peptide(dfs):
    # df_joined = reduce(lambda left,right: pd.merge(left, right), dfs)
    df_joined = pd.concat(dfs, axis=1)
    return df_joined


def read_and_join_files(mode, pep_colname,
                  samp_groups, int_file,
                  tax_file=None, func_file=None,
                  tax_colname=None, func_colname=None):

    # intensity
    int = read_intensity_table(int_file, samp_groups, pep_colname)

    # start df list
    dfs = [int]

    if mode == 'tax' or mode == 'taxfn':
        tax = read_taxonomy_table(tax_file, pep_colname, tax_colname)
        dfs.append(tax)

    if mode == 'fn' or mode == 'taxfn':
        func = read_function_table(func_file, pep_colname, func_colname)
        dfs.append(func)

    dfs_joined = join_on_peptide(dfs)
    return(dfs_joined)


class SampleGroups:
    """
    flatten sample names list
    :param sample_names: dictionary of lists, where the keys are the group names and
    the values are lists of column names within that group
    :return:
    """

    def __init__(self, sample_names):

        # top level dictionary
        self.sample_names = sample_names

        # flatten sample column names
        self.all_intcols = list(itertools.chain(*list(sample_names.values())))

        # for defining pandas column data types on read in
        self.dict_numeric_cols = {x: np.float64 for x in self.all_intcols}

        # number of conditions
        self.ngrps = len(sample_names)

        # name of experimental groups
        # sort alphabetically, so it's deterministic
        self.grp_names = sorted(list(sample_names.keys()))

        # when calculating means, column names for means
        self.mean_names = [grp + "_mean" for grp in self.grp_names]

