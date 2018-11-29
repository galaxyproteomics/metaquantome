import pandas as pd

from metaquantome.util.check_args import function_check, tax_check
from metaquantome.util.constants import MISSING_VALUES


def read_and_join_files(mode, pep_colname_int, pep_colname_func, pep_colname_tax, samp_grps, int_file, tax_file=None,
                        func_file=None, func_colname=None, tax_colname=None):
    """
    todo: doc
    :param pep_colname_func:
    :param pep_colname_tax:
    :param mode:
    :param pep_colname_int:
    :param samp_grps:
    :param int_file:
    :param tax_file:
    :param func_file:
    :param func_colname:
    :param tax_colname:
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

    dfs_joined = join_on_peptide(dfs)
    return dfs_joined


def read_intensity_table(file, samp_grps, pep_colname_int):
    """

    :param file:
    :param samp_grps:
    :param pep_colname_int:
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
    Peptides with no annotation are kept, and assigned 32644 (ncbi id for unassigned)
    :param data_dir:
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
    todo:doc
    :param file:
    :param pep_colname_func:
    :return:
    """
    df = pd.read_table(file, sep="\t", index_col=pep_colname_func,
                       na_values=MISSING_VALUES)
    df_new = df[[func_colname]].copy()
    # drop nas
    df_new.dropna(inplace=True, axis=0)
    return df_new


def read_nopep_table(file, mode, samp_grps, func_colname=None, tax_colname=None):
    """

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
    sub = list()
    if mode == 'f':
        sub = [func_colname]
    elif mode == 't':
        sub = [tax_colname]
    elif mode == 'ft':
        sub = [func_colname, tax_colname]

    df.dropna(how='all', subset=sub, inplace=True)

    # type_change = {col: object for col in sub}
    # df_new = df.astype(dtype=type_change)
    return df


def join_on_peptide(dfs):
    # todo: doc
    # join inner means that only peptides present in all dfs will be kept
    df_all = dfs.pop(0)
    while len(dfs) > 0:
        df_other = dfs.pop(0)
        df_all = df_all.join(df_other, how="inner")
    return df_all


def write_out_general(df, outfile, cols):
    # todo: doc
    df.to_csv(outfile,
              columns=cols,
              sep="\t",
              header=True,
              index=False,
              na_rep="NA")


def define_outfile_cols_expand(samp_grps, ontology, mode):
    # todo: doc
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
