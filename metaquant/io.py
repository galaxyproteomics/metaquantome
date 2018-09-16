import os

import pandas as pd

MISSING_VALUES = ["", "0", "NA", "NaN", "0.0"]
ONTOLOGIES = ['cog', 'go', 'ec']


def read_and_join_files(mode, pep_colname,
                        samp_groups, int_file,
                        tax_file=None, func_file=None,
                        func_colname=None,
                        tax_colname=None):

    # intensity
    int = read_intensity_table(int_file, samp_groups, pep_colname)

    # start df list
    dfs = [int]
    if mode == 'tax' or mode == 'taxfn':
        tax_check(tax_file, tax_colname)
        tax = read_taxonomy_table(tax_file, pep_colname, tax_colname)
        dfs.append(tax)
    if mode == 'fn' or mode == 'taxfn':
        function_check(func_file, func_colname)
        func = read_function_table(func_file, pep_colname, func_colname)
        dfs.append(func)

    dfs_joined = join_on_peptide(dfs)
    return dfs_joined


def read_intensity_table(file, samp_grps, pep_colname):
    # read in data
    df = pd.read_table(file, sep="\t", index_col=pep_colname,
                       dtype=samp_grps.dict_numeric_cols,
                       na_values=MISSING_VALUES,
                       low_memory=False)

    # drop rows where all intensities are NA
    df.dropna(axis=0, how="all", inplace=True)

    # change remaining missing intensities to 0, for arithmetic (changed back to NA for export)
    values = {x: 0 for x in samp_grps.all_intcols}
    df.fillna(values, inplace=True)

    return df


def read_taxonomy_table(file, pep_colname, tax_colname):
    """
    read taxonomy table, such as Unipept output.
    Peptides with no annotation are kept, and assigned 32644 (ncbi id for unassigned)
    :param data_dir:
    :param file: path to taxonomy file
    :param pep_colname: string, peptide sequence column name
    :param tax_colname: string, taxonomy identifier column name
    :return: a pandas dataframe where index is peptide sequence and the single column is the associated ncbi taxid
    """
    # always read as character
    df = pd.read_table(file, sep="\t", index_col=pep_colname,
                       na_values=MISSING_VALUES, dtype={tax_colname: object})

    # take only specified column
    df_tax = df.loc[:, [tax_colname]]
    return df_tax


def read_function_table(file, pep_colname, func_colname):
    """

    :param file:
    :param pep_colname:
    :return:
    """
    df = pd.read_table(file, sep="\t", index_col=pep_colname,
                       na_values=MISSING_VALUES)

    df_new = df[[func_colname]].copy()

    # drop nas
    df_new.dropna(inplace=True, axis=0)
    return(df_new)


def join_on_peptide(dfs):
    # join inner means that only peptides present in all dfs will be kept
    df_joined = pd.concat(dfs, axis=1, join='inner')
    return df_joined


def function_check(func_file, func_colname):
    if not func_file:
        raise IOError('Function tabular file not provided (--func_file)')
    if not os.path.exists(func_file):
        raise FileNotFoundError('func_file does not exist. Please check filename')
    if not func_colname:
        raise ValueError('func_colname=None. Please provide a function column name (--func_colname)')


def tax_check(tax_file, tax_colname):
    if not tax_file:
        raise IOError('Taxonomy tabular file not provided (--tax_file)')
    if not os.path.exists(tax_file):
        raise FileNotFoundError('func_file does not exist. Please check filename')
    if not tax_colname:
        raise ValueError('tax_colname=None. Please provide a taxonomy column name (--tax_colname)')