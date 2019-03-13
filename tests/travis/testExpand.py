import unittest
import numpy as np
import pandas as pd

from metaquantome.databases import GeneOntologyDb as godb
import metaquantome.modules.expand as expand
from metaquantome.classes.SampleGroups import SampleGroups
from metaquantome.util.testutils import testfile
from metaquantome.util.utils import TEST_DIR


class TestExpandUtils(unittest.TestCase):
    # test dataframe
    # 9604 has no peptides observed in either sample, but has enough from the other peptides
    # 9606 has one peptide observed in one sample and two in the other
    # 9605 has two peptides in each sample, but has 1 child and is a leaf
    # 9599 only has one peptide
    peptide = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    lca = [9604, 9604, 9605, 9605, 9606, 9606, 9599]
    samp1 = [1, 1, 1, 1, 1, 0, 2]
    samp2 = [0, 1, 1, 1, 1, 1, 2]
    samp3 = [1, 1, 0, 1, 1, 0, 2]
    samp_grps = SampleGroups('{"grp1": ["samp1", "samp2", "samp3"]}')
    test_df = pd.DataFrame({'lca': lca,
                            'samp1': samp1,
                            'samp2': samp2,
                            'samp3': samp3},
                           index=[peptide])

    def testCalcMeans(self):
        calk = expand.calc_means(self.test_df, self.samp_grps)['grp1_mean']
        expected_means = pd.Series(np.log2([2/3, 1, 2/3, 1, 1, 1/3, 2]),
                                   index=[self.peptide],
                                   name='grp1_mean')
        self.assertTrue(calk.equals(expected_means))

        df = pd.DataFrame({'s1_1': [4, 4],
                           's1_2': [2, 2],
                           's2_1': [5, 10],
                           's2_2': [7, 16]})
        samps = SampleGroups('{"s1": ["s1_1", "s1_2"], "s2": ["s2_1", "s2_2"]}')
        means = expand.calc_means(df, samps)
        self.assertTrue(means['s1_mean'].equals(pd.Series({0: np.log2(3.0), 1: np.log2(3.0)}, name="s1_mean")))


class TestFunctionalAnalysisExpand(unittest.TestCase):

    db = godb.GeneOntologyDb(TEST_DIR, slim_down=True)

    def testSingleInt(self):
        func=testfile('simple_func.tab')
        int=testfile('simple_int.tab')
        go_df = expand.expand('f', sinfo='{"s1": ["int"]}', int_file=int, pep_colname_int='peptide',
                              pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, func_file=func,
                              func_colname='go', ontology='go')
        self.assertEqual(go_df.loc["GO:0022610"]['int'], np.log2(200))
        self.assertEqual(go_df.loc["GO:0008152"]['int'], np.log2(100))

    def testMultipleInt(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        go_df = expand.expand('f', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=int, pep_colname_int='peptide',
                              pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, func_file=func,
                              func_colname='go', ontology='go')
        self.assertEqual(go_df.loc['GO:0008152']['int1'], np.log2(10))
        self.assertEqual(go_df.loc['GO:0022610']['int2'], np.log2(30))
        # missing values (zeros, nans, NA's, etc) are turned into NaN's
        self.assertTrue(np.isnan(go_df.loc['GO:0000003']['int3']))
        return go_df

    def testNopep(self):
        nopep=testfile('nopep.tab')
        go_df = expand.expand('f', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=None, pep_colname_int='peptide',
                              pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, func_colname='go',
                              ontology='go', nopep=True, nopep_file=nopep).sort_index(axis=1)
        self.assertEqual(go_df.loc['GO:0008152']['int1'], np.log2(10))
        self.assertEqual(go_df.loc['GO:0022610']['int2'], np.log2(30))
        # missing values (zeros, nans, NA's, etc) are turned into NaN's
        self.assertTrue(np.isnan(go_df.loc['GO:0000003']['int3']))
        # now, test that the results are the same as obtained through the peptide method
        df = self.testMultipleInt().sort_index(axis=1)
        self.assertTrue(df.equals(go_df))

    def testSlimDown(self):
        func=testfile('func_eggnog.tab')
        int=testfile('int_eggnog.tab')
        sinfo='{"NS": ["int737NS", "int852NS", "int867NS"], "WS": ["int737WS", "int852WS", "int867WS"]}'
        go_df = expand.expand('f', sinfo=sinfo, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
                              pep_colname_tax='peptide', data_dir=TEST_DIR, func_file=func, func_colname='go', ontology='go',
                              slim_down=True)
        # test that all go terms are in slim
        # load slim
        returned_gos = set(go_df['id'])
        # potential of unknown, so just drop that
        returned_gos.discard('unknown')
        self.assertTrue(returned_gos.issubset(self.db.goslim.keys()))

    def testCog(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        cog_df = expand.expand('f', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=int, pep_colname_int='peptide',
                               pep_colname_func='peptide', pep_colname_tax='peptide', func_file=func, func_colname='cog',
                               ontology='cog')
        self.assertEqual(cog_df.loc["C"]['s1_mean'], np.log2((10+20+70)/3))
        self.assertEqual(cog_df.loc["N"]['int2'], np.log2(30))

    def testSimpleEc(self):
        func=testfile('simple_ec.tab')
        int=testfile('simple_int.tab')
        ec_df = expand.expand('f', sinfo='{"s1": ["int"]}', int_file=int, pep_colname_int='peptide',
                              pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, func_file=func,
                              func_colname='ec', ontology='ec')
        self.assertEqual(ec_df.loc["3.4.11.-"]['int'], np.log2(100))
        self.assertEqual(ec_df.loc["3.4.-.-"]['int'], np.log2(300))

    def testMultipleEc(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        ec_df = expand.expand('f', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=int, pep_colname_int='peptide',
                              pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, func_file=func,
                              func_colname='ec', ontology='ec')
        self.assertEqual(ec_df.loc['3.4.-.-']['int1'], np.log2(50))
        self.assertEqual(ec_df.loc['1.2.-.-']['int2'], np.log2(50))
        # missing values (zeros, nans, NA's, etc) are turned into NaN's
        self.assertTrue(np.isnan(ec_df.loc['1.2.-.-']['int3']))


class TestTaxonomyAnalysisExpand(unittest.TestCase):
    def testSingleBasic(self):
        tax = testfile('simple_tax.tab')
        int = testfile('simple_int.tab')
        tax_df = expand.expand('t', sinfo='{"s1": ["int"]}', int_file=int, pep_colname_int='peptide',
                               pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, tax_file=tax,
                               tax_colname='lca')
        self.assertEqual(tax_df.query("taxon_name == 'Helicobacter pylori'")['int'].values, np.log2(100))

    def testWrite(self):
        tax = testfile('simple_tax.tab')
        int = testfile('simple_int.tab')
        out = testfile('taxonomy_write_simple.tab')
        df = expand.expand(mode='t', sinfo='{"samp1": ["int"]}', int_file=int, pep_colname_int='peptide',
                           pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, outfile=out, tax_file=tax,
                           tax_colname='lca')
        written = pd.read_table(out)
        self.assertAlmostEqual(written.query("taxon_name == 'Clostridioides difficile'")['samp1_mean'].values[0], np.log2(200))

    def testMultCols(self):
        tax=testfile('multiple_tax.tab')
        int=testfile('multiple_int.tab')
        tax_df = expand.expand('t', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=int, pep_colname_int='peptide',
                               pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, tax_file=tax,
                               tax_colname='lca')
        self.assertEqual(tax_df.query("rank == 'phylum' and taxon_name == 'Proteobacteria'")['int3'].values[0], np.log2(70))

    def testNopep(self):
        nopep=testfile('nopep.tab')
        tax_df = expand.expand('t', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=None, pep_colname_int='peptide',
                               pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, tax_colname='lca',
                               nopep=True, nopep_file=nopep)
        self.assertEqual(tax_df.query("rank == 'phylum' and taxon_name == 'Proteobacteria'")['int3'].values[0],
                         np.log2(70))

    def testParentIntensityHigher(self):
        """
        make sure that parents always have higher intensity than children
        """
        tax=testfile('test_root_sum_uni.tab')
        int=testfile('test_root_sum_int.tab')
        tax_df = expand.expand('t', sinfo='{"A": ["int"]}', int_file=int, pep_colname_int='peptide',
                               pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, tax_file=tax,
                               tax_colname='taxon_id')
        # filter to phylum and below
        tax_df_filt = tax_df[(tax_df["rank"] != 'no rank') & (tax_df["rank"] != 'superkingdom')]
        # firmicutes phylum should be highest
        ints = tax_df_filt['int']
        self.assertEqual(ints.max(), ints[1239])
        # strep genus intensity should be greater than or equal to that of strep species
        self.assertGreaterEqual(ints[1301], ints[1302])
        self.assertGreaterEqual(ints[1301], ints[1305])


class TestFunctionTaxonomyAnalysis(unittest.TestCase):
    def testDifferentNames(self):
        tax = testfile('ft_tax.tab')
        func = testfile('ft_func.tab')
        int = testfile('ft_int.tab')
        ft = expand.expand('ft', sinfo='{"A": ["int"]}', int_file=int, pep_colname_int='Sequence',
                           pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, tax_file=tax,
                           tax_colname='lca', func_file=func, func_colname="go")
        self.assertIn("A_mean", list(ft))


if __name__=='__main__':
    unittest.main()
