import pandas as pd
import numpy as np


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

    def __init__(self, file, intcols):
        numeric_cols = {x: np.float64 for x in intcols}

        # read in data
        df = pd.read_table(file, sep="\t", index_col="peptide", dtype = numeric_cols)

        # drop columns where all are NA
        df.dropna(axis=1, how="all", inplace = True)

        # change missing intensities to 0
        df[intcols] = df[intcols].fillna(0)

        # change missing taxa to 'unknown', because nan compares false to itself
        df.fillna('unknown', inplace=True)

        self.input = df

        # determine which tax. ranks are in data
        self.local_tax =  set(self.tax).intersection(set(self.input))

        # for later
        self.intcols = intcols

    def rank_contents(self, ranks):
        return set(self.input[ranks].values)

    def rel_abundance_rank(self, rank):
        summed_abund = self.input.groupby(by=rank)[self.intcols].sum(axis = 0)
        rel_abundance = summed_abund / summed_abund.sum(axis = 0)
        rel_abundance['rank'] = rank
        rel_abundance['member'] = rel_abundance.index
        if len(set(self.intcols)) > 1:
            new_columns = {x : 'rsra_' + x for x in self.intcols}
            rel_abundance.rename(mapper = new_columns, axis = 'columns', inplace = True)
        else:
            rel_abundance.rename(mapper = {self.intcols[0]: 'rsra'}, axis='columns', inplace = True)
        return rel_abundance

    def rel_abundance_all_ranks(self):
        all = pd.concat([self.rel_abundance_rank(x) for x in self.local_tax])
        return all

    def write_out(self, outfile):
        all = self.rel_abundance_all_ranks()
        all.to_csv(outfile, sep="\t", index=False)

# TODO
# CLI
