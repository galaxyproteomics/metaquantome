import pandas as pd
import numpy as np
import scipy.stats as sps
import statsmodels.sandbox.stats.multicomp as mc
from src import common


def taxonomy_analysis(file,
                      sample1_colnames,
                      sample2_colnames=None,
                      test=False,
                      threshold=2,
                      paired=False,
                      outfile=None):

    # define intensity columns
    all_intcols, dict_numeric_cols = common.define_intensity_columns(sample1_colnames, sample2_colnames)

    # read in data
    df = common.read_data_table(file, dict_numeric_cols, all_intcols)

    # determine which taxa are in user-provided dataset
    user_tax = set(tax).intersection(set(df))

    # add up through ranks
    norm_intensity_all_ranks = pd.concat([rel_abundance_rank(df, x, all_intcols) for x in user_tax])

    # test
    if test:
        results = test_norm_intensity(norm_intensity_all_ranks, sample1_colnames, sample2_colnames, threshold, paired)
    else:
        results = norm_intensity_all_ranks

    if outfile:
        write_out(results, outfile)
    else:
        return results


tax = ["superkingdom",
       "kingdom",
       "subkingdom",
       "superphylum",
       "phylum",
       "subphylum",
       "superclass",
       "class",
       "subclass"
       "infraclass",
       "superorder",
       "order",
       "suborder",
       "infraorder",
       "parvorder",
       "superfamily",
       "family",
       "subfamily",
       "tribe",
       "subtribe",
       "genus",
       "subgenus",
       "species_group",
       "species_subgroup",
       "species",
       "subspecies",
       "varietas",
       "forma"]


def rel_abundance_rank(df, rank, all_intcols):
    summed_abund = df.groupby(by=rank)[all_intcols].sum(axis=0)
    rel_abundance = summed_abund / summed_abund.sum(axis=0)
    rel_abundance['rank'] = rank
    rel_abundance['member'] = rel_abundance.index
    return rel_abundance


def test_norm_intensity(df, sample1_colnames, sample2_colnames, threshold, paired):
    # filter to a minimum number of observed intensities per sample
    numNotNA1 = (df[sample1_colnames] > 0).apply(sum, axis=1) >= threshold
    numNotNA2 = (df[sample2_colnames] > 0).apply(sum, axis=1) >= threshold
    keep = numNotNA1 & numNotNA2
    filtered_df = df.loc[keep].copy()

    # add mean, sd; create dataframe
    test_results = filtered_df.apply(lambda x: sps.stats.ttest_ind(x[sample1_colnames],
                                                                   x[sample2_colnames],
                                                                   equal_var=paired).pvalue, axis=1)
    filtered_df['log2ratio_grp2_over_grp1'] = \
        filtered_df.apply(
            lambda x: np.log2(np.mean(x[sample2_colnames]) / np.mean(x[sample1_colnames])),
        axis=1)
    filtered_df['p'] = test_results
    filtered_df['corrected_p'] = mc.fdrcorrection0(test_results, method='indep')[1]
    return filtered_df


def write_out(df, outfile):
    df.to_csv(outfile, sep="\t", index=False)
