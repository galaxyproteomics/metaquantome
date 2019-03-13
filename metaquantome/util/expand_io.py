import pandas as pd

from metaquantome.util.check_args import function_check, tax_check
from metaquantome.util.utils import MISSING_VALUES


def read_and_join_files(mode, pep_colname_int, pep_colname_func, pep_colname_tax, samp_grps, int_file, tax_file=None,
                        func_file=None, func_colname=None, tax_colname=None):
    """
    Reads in intensity, function, and/or taxonomy files, and joins all on the peptide column.

    :param pep_colname_func: name of the peptide column in the function file
    :param pep_colname_tax: name of the peptide column in the taxonomy file
    :param pep_colname_int: name of the peptide column in the intensity file
    :param mode: analysis mode - either 'f', 't', or 'ft'
    :param samp_grps: SampleGroups() object
    :param int_file: path to intensity file
    :param tax_file: path to taxonomy file. required for 't' and 'ft' modes
    :param func_file: path to function file. required for 'f' and 'ft' modes
    :param func_colname: column name of functional annotation in function file
    :param tax_colname: column name of taxonomic annotation in taxonomy file
    :return: joined dataframe; missing intensities as 0.
    """

    # intensity
    int = read_intensity_table(int_file, samp_grps, pep_colname_int)

    # start df list
    dfs = [int]
    if mode == 't' or mode == 'ft':
        tax_check(tax_file, tax_colname)
        tax = read_taxonomy_table(tax_file, pep_colname_tax, tax_colname)
        dfs.append(tax)
    if mode == 'f' or mode == 'ft':
        function_check(func_file, func_colname)
        func = read_function_table(func_file, pep_colname_func, func_colname)
        dfs.append(func)
    # join all
    dfs_joined = join_on_peptide(dfs)
    dfs_joined.index.name = 'peptide'
    return dfs_joined


def read_intensity_table(file, samp_grps, pep_colname_int):
    """
    read the file containing peptide intensities to a pandas dataframe.

    :param file: path to intensity file. must be tab-separated
    :param samp_grps: SampleGroups object
    :param pep_colname_int: name of peptide column in intensity table
    :return: intensity table; missing values as 0
    """
    # read in data
    df = pd.read_table(file, sep="\t", index_col=pep_colname_int,
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


def read_taxonomy_table(file, pep_colname_tax, tax_colname):
    """
    read taxonomy table, such as Unipept output.
    Peptides with no annotation are dropped.

    :param file: path to taxonomy file
    :param pep_colname_tax: string, peptide sequence column name
    :param tax_colname: string, taxonomy identifier column name
    :return: a pandas dataframe where index is peptide sequence and the single column is the associated ncbi taxid
    """
    # always read as character
    df = pd.read_table(file, sep="\t", index_col=pep_colname_tax,
                       na_values=MISSING_VALUES, dtype={tax_colname: object})
    # take only specified column
    df_tax = df.loc[:, [tax_colname]]

    # drop nas
    df_tax.dropna(inplace=True, axis=0)
    return df_tax


def read_function_table(file, pep_colname_func, func_colname):
    """
    read functional annotation table to Pandas dataframe. Peptides
    with no annotation are dropped.

    :param file: path to tab-separated function file
    :param pep_colname_func: name of peptide column in function column
    :param func_colname: name of functional annotation column in function table
    :return: pandas dataframe where index is peptide sequence and single column is associated functional annotation.
    """
    df = pd.read_table(file, sep="\t", index_col=pep_colname_func,
                       na_values=MISSING_VALUES)
    df_new = df[[func_colname]].copy()
    # drop nas
    df_new.dropna(inplace=True, axis=0)
    return df_new


def join_on_peptide(dfs):
    """
    Inner join a list of dataframes on the index.

    :param dfs: list of pandas dataframes
    :return: joined dataframe.
    """
    # join inner means that only peptides present in all dfs will be kept
    df_all = dfs.pop(0)
    while len(dfs) > 0:
        df_other = dfs.pop(0)
        df_all = df_all.join(df_other, how="inner")
    return df_all


def read_nopep_table(file, mode, samp_grps, func_colname=None, tax_colname=None):
    """
    Read in a pre-joined table (rather than 3 separate tables)

    :param file: file with intensity and functional or taxonomic terms
    :param mode: f, tax, or ft
    :param samp_grps: SampleGroups() object
    :param func_colname: name of column with functional terms
    :param tax_colname: name of column with taxonomic annotations
    :return: dataframe, missing values as 0
    """
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

    # drop rows where function is missing (for mode 'f'), taxonomy is missing (mode 't'),
    # or both function and taxonomy are missing (mode 'ft')
    sub = list()
    if mode == 'f':
        sub = [func_colname]
    elif mode == 't':
        sub = [tax_colname]
    elif mode == 'ft':
        sub = [func_colname, tax_colname]
    df.dropna(how='all', subset=sub, inplace=True)
    return df


def write_out_general(df, outfile, cols):
    """
    Write a pandas dataframe as a tab-separated file.
    Keeps header, does not write index; missing
    values are represented as NA

    :param df: dataframe
    :param outfile: path to output file
    :param cols: columns to be written, in desired order
    :return: None
    """
    df.to_csv(outfile,
              columns=cols,
              sep="\t",
              header=True,
              index=False,
              na_rep="NA")


def define_outfile_cols_expand(samp_grps, ontology, mode):
    """
    define columns for writing the expand output file

    :param samp_grps: SampleGroups object
    :param ontology: functional ontology. only required for 'f' or 'ft' modes
    :param mode: f, t, or ft
    :return: a list of relevant columns in the correct order
    """
    int_cols = []
    int_cols += samp_grps.mean_names + samp_grps.all_intcols
    node_cols = []
    if ontology != "cog":
        node_cols += samp_grps.n_peptide_names_flat
        # ft doesn't have samp_children
        if mode != 'ft':
            node_cols += samp_grps.samp_children_names_flat
    quant_cols = int_cols + node_cols
    if mode == 'f':
        if ontology == 'go':
            cols = ['id', 'name', 'namespace'] + quant_cols
        elif ontology == 'cog':
            cols = ['id', 'description'] + quant_cols
        elif ontology == 'ec':
            cols = ['id', 'description'] + quant_cols
        else:
            raise ValueError("Invalid ontology. Expected one of: %s" % ['go', 'cog', 'ec'])
    elif mode == 't':
        cols = ['id', 'taxon_name', 'rank'] + quant_cols
    elif mode == 'ft':
        cols = ['go_id', 'name', 'namespace', 'tax_id', 'taxon_name', 'rank'] + quant_cols
    else:
        raise ValueError("Invalid mode. Expected one of: %s" % ['f', 't', 'ft'])
    return cols
