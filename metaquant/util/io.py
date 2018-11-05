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
    # only intcols (in case table has extra cols)
    int_df = df.loc[:, samp_grps.all_intcols]

    # drop rows where all intensities are NA
    int_df.dropna(axis=0, how="all", inplace=True)
    # change remaining missing intensities to 0, for arithmetic (changed back to NA for export)
    values = {x: 0 for x in samp_grps.all_intcols}
    int_df.fillna(values, inplace=True)
    return int_df


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

    # drop nas
    df_tax.dropna(inplace=True, axis=0)
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
    return df_new


def read_nopep_table(file, mode, samp_grps, func_colname=None, tax_colname=None):
    newdict = samp_grps.dict_numeric_cols.copy()
    newdict[func_colname] = object
    newdict[tax_colname] = object
    df = pd.read_table(file, sep="\t",
                       dtype=newdict,
                       na_values=MISSING_VALUES,
                       low_memory=False)
    # change remaining missing intensities to 0, for arithmetic (changed back to NA for export)
    values = {x: 0 for x in samp_grps.all_intcols}
    df.fillna(values, inplace=True)
    sub = list()
    if mode == 'fn':
        sub = [func_colname]
    elif mode == 'tax':
        sub = [tax_colname]
    elif mode == 'taxfn':
        sub = [func_colname, tax_colname]

    df.dropna(how='all', subset=sub, inplace=True)

    # type_change = {col: object for col in sub}
    # df_new = df.astype(dtype=type_change)
    return df



def join_on_peptide(dfs):
    # join inner means that only peptides present in all dfs will be kept
    df_all = dfs.pop(0)
    while len(dfs) > 0:
        df_other = dfs.pop(0)
        df_all = df_all.join(df_other, how="inner")
    return df_all


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


def define_outfile_cols_expand(samp_grps, ontology, mode):
    int_cols = []
    int_cols += samp_grps.mean_names + samp_grps.all_intcols
    if mode == 'fn':
        if ontology == 'go':
            cols = ['id', 'name', 'namespace'] + int_cols
        elif ontology == 'cog':
            cols = ['id', 'description'] + int_cols
        elif ontology == 'ec':
            cols = ['id', 'description'] + int_cols
        else:
            raise ValueError("Invalid ontology. Expected one of: %s" % ['go', 'cog', 'ec'])
    elif mode == 'tax':
        cols = ['id', 'taxon_name', 'rank'] + int_cols
    elif mode == 'taxfn':
        cols = [ontology, 'cog_descript', 'taxon_name', 'rank'] + int_cols
    else:
        raise ValueError("Invalid mode. Expected one of: %s" % ['fun', 'tax', 'taxfn'])
    return cols
