import pandas as pd

from metaquantome.util.utils import MISSING_VALUES, P_COLNAME, P_CORR_COLNAME
import metaquantome.util.expand_io as expand_io


def read_expanded_table(file, samp_grps):
    """
    read the output of metaquantome.expand

    :param file: path to file
    :param samp_grps: SampleGroups object
    :return: dataframe, with missing values represented as 0
    """
    df = pd.read_table(file, sep="\t",
                       dtype=samp_grps.dict_numeric_cols_expanded,
                       na_values=MISSING_VALUES,
                       low_memory=False)
    df.fillna(0, inplace=True)
    return df


def write_stat(df, outfile, samp_grps, ontology, mode):
    """
    write the output of stat

    :param df: data frame
    :param outfile: path to output file
    :param samp_grps: SampleGroups object
    :param ontology: functional ontology (for 'f' or 'ft' modes)
    :param mode: f, t, or ft
    :return: None
    """
    cols = expand_io.define_outfile_cols_expand(samp_grps, ontology, mode) +\
           [samp_grps.fc_name, P_COLNAME, P_CORR_COLNAME]
    expand_io.write_out_general(df, outfile=outfile, cols=cols)
