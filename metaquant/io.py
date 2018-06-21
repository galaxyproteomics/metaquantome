import pandas as pd
import itertools
from metaquant import phylo_tree
import re
import numpy as np
import json
import os


MISSING_VALUES = ["", "0", "NA", "NaN", "0.0"]
ONTOLOGIES = ['cog', 'go']


def read_intensity_table(file, samp_grps, pep_colname):
    # read in data
    df = pd.read_table(file, sep="\t", index_col=pep_colname,
                       dtype=samp_grps.dict_numeric_cols,
                       na_values=MISSING_VALUES,
                       low_memory=False)

    # drop columns where all are NA...which sometimes happens for whatever reason
    df.dropna(axis=1, how="all", inplace=True)

    # change missing intensities to 0, for arithmetic (changed back to NA for export)
    values = {x: 0 for x in samp_grps.all_intcols}
    df.fillna(values, inplace=True)

    return df


def read_taxonomy_table(file, pep_colname, tax_colname):
    """
    read taxonomy table, such as Unipept output.
    Peptides with no annotation are kept, and assigned 32644 (ncbi id for unassigned)
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

    # check for numeric characters, which indicates taxid
    # if is name, convert to taxid
    # keep as character until querying ncbi database
    if sniff_tax_names(df_tax, tax_colname):
        ncbi = phylo_tree.load_ncbi()
        df_tax[tax_colname] = phylo_tree.convert_name_to_taxid(df_tax[tax_colname], ncbi)

    return df_tax


def sniff_tax_names(df, tax_colname):
    pattern = re.compile(r'[0-9]')  # little bit faster to compile
    is_numeric = df[tax_colname].str.contains(pattern)
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
    # join inner means that only peptides present in all dfs will be kept
    df_joined = pd.concat(dfs, axis=1, join='inner')
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
        # same order as grp names
        self.mean_names = [grp + "_mean" for grp in self.grp_names]


def read_samp_info(sinfo):
    # check if sinfo is a file name
    if os.path.exists(sinfo):
        samp_names = dict()
        with open(sinfo, 'r') as f:
            # throw away header
            f.readline()

            # read groups one by one
            for line in f:
                split = line.split('\t')
                samp_names[split[0]] = [elem.strip() for elem in split[1].split(',')]
        return samp_names
    else:
        # check if sinfo is a json format
        is_json, json_obj = to_json(sinfo)

        if is_json:
            return json_obj
        else:
            raise ValueError('--samps is not a text file or in proper json format. please check again!')


# thanks to https://stackoverflow.com/questions/5508509/how-do-i-check-if-a-string-is-valid-json-in-python
def to_json(obj):
    try:
        json_object = json.loads(obj)
        return True, json_object
    except ValueError:
        return False, obj