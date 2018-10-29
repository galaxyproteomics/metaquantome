import unittest
import numpy as np
import os

from metaquant.databases import GeneOntologyDb
from metaquant.analysis.expand import expand
from metaquant.analysis.test import test
from tests.testutils import testfile, TTEST_SINFO
from metaquant.util.utils import DATA_DIR, GO_SUBDIR


class TestFunctionalAnalysisExpand(unittest.TestCase):
    def testSingleInt(self):
        func=testfile('simple_func.tab')
        int=testfile('simple_int.tab')
        go_df = expand('fn', samps='{"s1": ["int"]}', int_file=int, pep_colname='peptide', func_file=func,
                       func_colname='go', ontology='go')
        self.assertEqual(go_df.loc["GO:0022610"]['int'], np.log2(200))
        self.assertEqual(go_df.loc["GO:0008152"]['int'], np.log2(100))

    def testMultipleInt(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        go_df = expand('fn', samps='{"s1": ["int1", "int2", "int3"]}', int_file=int, func_colname='go', func_file=func,
                       ontology='go')
        self.assertEqual(go_df.loc['GO:0008152']['int1'], np.log2(10))
        self.assertEqual(go_df.loc['GO:0022610']['int2'], np.log2(30))
        # missing values (zeros, nans, NA's, etc) are turned into NaN's
        self.assertTrue(np.isnan(go_df.loc['GO:0000003']['int3']))

    def testSlimDown(self):
        func=testfile('func_eggnog.tab')
        int=testfile('int_eggnog.tab')
        sinfo='{"NS": ["int737NS", "int852NS", "int867NS"], "WS": ["int737WS", "int852WS", "int867WS"]}'
        go_df = expand('fn', samps=sinfo, int_file=int, func_colname='go', func_file=func, ontology='go',
                       slim_down=True)
        # test that all go terms are in slim
        # load slim
        go = GeneOntologyDb.GeneOntologyDb(data_dir=os.path.join(DATA_DIR, GO_SUBDIR), slim_down=True, overwrite=False)
        returned_gos = set(go_df['id'])
        # potential of unknown, so just drop that
        returned_gos.discard('unknown')
        self.assertTrue(returned_gos.issubset(go.goslim.keys()))

    def testCog(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        cog_df = expand('fn', samps='{"s1": ["int1", "int2", "int3"]}', int_file=int, func_colname='cog',
                        func_file=func, ontology='cog')
        self.assertEqual(cog_df.loc["C"]['s1_mean'], np.log2((10+20+70)/3))
        self.assertEqual(cog_df.loc["N"]['int2'], np.log2(30))

    def testSimpleEc(self):
        func=testfile('simple_ec.tab')
        int=testfile('simple_int.tab')
        ec_df = expand('fn', samps='{"s1": ["int"]}', int_file=int, pep_colname='peptide', func_colname='ec',
                       func_file=func, ontology='ec')
        self.assertEqual(ec_df.loc["3.4.11.-"]['int'], np.log2(100))
        self.assertEqual(ec_df.loc["3.4.-.-"]['int'], np.log2(300))

    def testMultipleEc(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        ec_df = expand('fn', samps='{"s1": ["int1", "int2", "int3"]}', int_file=int, func_colname='ec', func_file=func,
                                 ontology='ec')
        self.assertEqual(ec_df.loc['3.4.-.-']['int1'], np.log2(50))
        self.assertEqual(ec_df.loc['1.2.-.-']['int2'], np.log2(50))
        # missing values (zeros, nans, NA's, etc) are turned into NaN's
        self.assertTrue(np.isnan(ec_df.loc['1.2.-.-']['int3']))


class TestFunctionalAnalysisTest(unittest.TestCase):
    def testDA(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        df_expd = expand('fn', samps=TTEST_SINFO, int_file=int, func_colname='go', func_file=func, ontology='go')
        df_tst = test(df_expd, samps=TTEST_SINFO, paired=False, parametric=True)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(df_tst['p']['GO:0008152'] > 0.05)
        self.assertTrue(df_tst['p'][['GO:0022610','GO:0000003','GO:0032505']].le(0.05).all())

    def testCogTTest(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        cog_df = expand('fn', samps=TTEST_SINFO, int_file=int, func_colname='cog', func_file=func, ontology='cog')
        cog_tst = test(cog_df, samps=TTEST_SINFO, paired=False, parametric=True)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(cog_df['p']['C'] > 0.05)
        self.assertTrue(cog_df['p'][['N','D']].le(0.05).all())

    def testDiffAbundEc(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        ec_df = expand('fn', samps=TTEST_SINFO, int_file=int, func_colname='ec', func_file=func, ontology='ec')
        ec_tst = test(ec_df, samps=TTEST_SINFO, paired=False, parametric=True)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(ec_df['p']['3.4.11.-'] > 0.05)
        self.assertTrue(ec_df['p'][['3.4.21.70','1.2.-.-']].le(0.05).all())


if __name__ == '__main__':
    unittest.main()
