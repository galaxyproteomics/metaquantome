import pandas as pd

from metaquantome.util.constants import MISSING_VALUES, P_COLNAME, P_CORR_COLNAME
import metaquantome.util.expand_io as expand_io


def read_expanded_table(file, samp_grps):
    df = pd.read_table(file, sep="\t",
                       dtype=samp_grps.dict_numeric_cols_expanded,
                       na_values=MISSING_VALUES,
                       low_memory=False)
    return df


def write_test(df, outfile, samp_grps, ontology, mode):
    # todo: doc
    cols = expand_io.define_outfile_cols_expand(samp_grps, ontology, mode) +\
           [samp_grps.fc_name, P_COLNAME, P_CORR_COLNAME]
    expand_io.write_out_general(df, outfile=outfile, cols=cols)
