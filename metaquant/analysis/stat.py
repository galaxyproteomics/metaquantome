import pandas as pd
import sys

from metaquant.util.io import MISSING_VALUES, define_outfile_cols_expand
from metaquant.SampleGroups import SampleGroups
from metaquant.util import stats
from metaquant.analysis import separation


def stat(infile, outfile, samps, mode, ontology, statmode='t', paired=0, parametric=0):
    df = read_expanded(infile)
    # define sample groups
    samp_grps = SampleGroups(samps)
    if statmode == 't':
        if samp_grps.ngrps != 2:
            ValueError('testing is only available for 2 experimental conditions.')
        # run test
        df_test = stats.test_norm_intensity(infile, samp_grps, paired, parametric)
        # write outfile
        write_test(df=df_test, samp_grps=samp_grps, ontology=ontology, mode=mode, outfile=outfile)
        # return
        return df_test
    elif statmode == 'sep':
        # todo: WIP
        df.fillna(0, inplace=True)
        # run sep
        sep = separation.separation_of_n_clusts(df, samp_grps)
        # print sep to stdout
        with open(outfile, 'w') as of:
            header = 'mode: ' + mode
            if mode != 'tax':
                header += '; ontology: ' + ontology
            header += '\n'
            of.write(header)
            of.write(str(sep))
        return sep
    else:
        ValueError('incorrect statmode. Should be "t" or "sep"')


def read_expanded(file):
    df = pd.read_table(file, sep="\t",
                       na_values=MISSING_VALUES)
    return df


def write_test(df, outfile, samp_grps, ontology, mode):
    cols = define_outfile_cols_expand(samp_grps, ontology, mode) +\
        [samp_grps.fc_name, stats.P_COLNAME, stats.P_CORR_COLNAME]
    df.to_csv(outfile, columns=cols, sep="\t", header=True, index=False, na_rep="NA")
