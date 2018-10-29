import pandas as pd

from metaquant.util.io import MISSING_VALUES, define_outfile_cols_expand
from metaquant.SampleGroups import SampleGroups
from metaquant.util.stats import test_norm_intensity


def test(df, samps, paired=0, parametric=0):

    # define sample groups
    samp_grps = SampleGroups(samps)

    if samp_grps.ngrps != 2:
        ValueError('testing is only available for 2 experimental conditions.')

    # run test
    df_test = test_norm_intensity(df, samp_grps, paired, parametric)

    # return
    return df_test


def read_expanded(file):
    df = pd.read_table(file, sep="\t",
                       na_values=MISSING_VALUES)
    return df