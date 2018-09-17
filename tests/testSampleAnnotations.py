import unittest
import pandas as pd
from metaquant.databases.NCBITaxonomyDb import NCBITaxonomyDb
from metaquant.util.utils import define_ontology_data_dir
from metaquant.SampleAnnotations import SampleAnnotations
from metaquant.SampleGroups import SampleGroups
import tests.testutils  # runs pandas options
import numpy as np


class TestSampleAnnotations(unittest.TestCase):
    ddir = define_ontology_data_dir('taxonomy')
    ncbi = NCBITaxonomyDb(ddir)

    # todo: add test init

    def testAddSample(self):
        peptide = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
        # 9604 has no peptides observed in either sample, but has enough from the other peptides
        # 9606 has one peptide observed in one sample and two in the other
        # 9605 has two peptides in each sample, but has 1 child and is a leaf
        # 9599 only has one peptide
        lca = [9604, 9604, 9605, 9605, 9606, 9606, 9599]
        samp1 = [1, 1, 1, 1, 1, 0, 1]
        samp2 = [1, 1, 1, 1, 1, 1, 1]
        samp_grps = SampleGroups('{"grp1": ["samp1", "samp2"]}')
        test_df = pd.DataFrame({'lca': lca,
                                'samp1': samp1,
                                'samp2': samp2},
                               index=[peptide])
        sa = SampleAnnotations(db=self.ncbi)
        sa.add_samples_from_df(test_df, annot_colname='lca', samp_grps=samp_grps,
                               min_peptides=2, min_children_non_leaf=2)
        df = sa.to_dataframe().sort_index(axis=1)
        exp_df = pd.DataFrame({'samp1': [6, np.nan],
                               'samp2': [7, 2]},
                              index=[9604, 9606]).sort_index(axis=1)
        self.assertTrue(df.equals(exp_df))


if __name__ == '__main__':
    unittest.main()