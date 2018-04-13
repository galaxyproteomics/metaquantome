import pandas as pd
import numpy as np
import scipy.stats as sps
import statsmodels.sandbox.stats.multicomp as mc

class TreeOfLife:
    """contains taxonomic ranks and intensities of associated members"""
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

    def __init__(self, file, grp1_intcols, grp2_intcols=None):

        self.grp1 = grp1_intcols
        self.grp2 = grp2_intcols

        if self.grp2:
            self.all_intcols = self.grp1 + self.grp2
        else:
            self.all_intcols = grp1_intcols

        numeric_cols = {x: np.float64 for x in self.all_intcols}

        # read in data
        df = pd.read_table(file, sep="\t", index_col="peptide", dtype = numeric_cols)

        # drop columns where all are NA
        df.dropna(axis=1, how="all", inplace = True)

        # change missing intensities to 0
        df[self.all_intcols] = df[self.all_intcols].fillna(0)

        # change missing taxa to 'unknown', because nan compares false to itself
        df.fillna('unknown', inplace=True)

        self.input = df

        # determine which tax. ranks are in data
        self.local_tax = set(self.tax).intersection(set(self.input))

    def rank_contents(self, ranks):
        return set(self.input[ranks].values)

    def rel_abundance_rank(self, rank):
        summed_abund = self.input.groupby(by=rank)[self.all_intcols].sum(axis = 0)
        rel_abundance = summed_abund / summed_abund.sum(axis = 0)
        rel_abundance['rank'] = rank
        rel_abundance['member'] = rel_abundance.index
        # if len(set(self.all_intcols)) > 1:
        #     new_columns = {x : 'rsra_' + x for x in self.all_intcols}
        #     rel_abundance.rename(mapper = new_columns, axis = 'columns', inplace = True)
        # else:
        #     rel_abundance.rename(mapper = {self.all_intcols[0]: 'rsra'}, axis='columns', inplace = True)
        return rel_abundance

    def rel_abundance_all_ranks(self):
        all = pd.concat([self.rel_abundance_rank(x) for x in self.local_tax])
        return all

    def test_norm_intensity(self):
        rel_abund = self.rel_abundance_all_ranks()
        # add mean, sd; create dataframe
        test_results = rel_abund.apply(lambda x: sps.stats.ttest_ind(x[self.grp1],
                                                                     x[self.grp2],
                                                                     equal_var = False).pvalue, axis = 1)
        corrected_p = mc.fdrcorrection0(test_results, method='indep')[1]
        return corrected_p

    def write_out(self, outfile):
        all = self.rel_abundance_all_ranks()
        all.to_csv(outfile, sep="\t", index=False)

# TODO
# add t.tests, fdr correction
