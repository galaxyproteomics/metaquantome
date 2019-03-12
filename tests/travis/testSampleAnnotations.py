import unittest
import pandas as pd

from metaquantome.databases.NCBITaxonomyDb import NCBITaxonomyDb
from metaquantome.classes.SampleAnnotations import SampleAnnotations
from metaquantome.classes.SampleGroups import SampleGroups
from metaquantome.util.utils import TEST_DIR


class TestSampleAnnotations(unittest.TestCase):
    ncbi = NCBITaxonomyDb(TEST_DIR)

    # test dataframe
    # 9604 has no peptides observed in either sample, but has enough from the other peptides
    # 9606 has one peptide observed in one sample and two in the other
    # 9605 has two peptides in each sample, but has 1 child and is a leaf
    # 9599 only has one peptide
    peptide = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    lca = [9604, 9604, 9605, 9605, 9606, 9606, 9599]
    samp1 = [1, 1, 1, 1, 1, 0, 1]
    samp2 = [1, 1, 1, 1, 1, 1, 1]
    samp_grps = SampleGroups('{"grp1": ["samp1", "samp2"]}')
    test_df = pd.DataFrame({'lca': lca,
                            'samp1': samp1,
                            'samp2': samp2},
                           index=[peptide])

    def testInit(self):
        init = SampleAnnotations(self.ncbi)
        self.assertEqual(self.ncbi, init.db)
        self.assertEqual(set(), init.hierarchies)

    def testAddSample(self):
        sa = SampleAnnotations(db=self.ncbi)
        sa.add_samples_from_df(self.test_df, annot_colname='lca',
                               samp_grps=self.samp_grps)
        df = sa.to_dataframe().sort_index(axis=1)
        self.assertEqual(df.loc[9604, 'samp1'], 6)
        self.assertEqual(df.loc[9599, 'samp1'], 1)

    def testCreateSampleSet(self):
        sampset = SampleAnnotations.create_sample_set(self.test_df, 'lca')
        self.assertSetEqual(sampset, set(self.lca))

    def testFilterToSampObserved(self):
        ftso = SampleAnnotations.filter_to_samp_observed(self.test_df, 'samp1')
        expected = self.test_df[~self.test_df.index.isin(['F'])]
        self.assertTrue(expected.equals(ftso))


if __name__ == '__main__':
    unittest.main()
