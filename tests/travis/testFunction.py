import unittest
import numpy as np


from metaquantome.databases import GeneOntologyDb as godb
from metaquantome.analysis.expand import expand
from metaquantome.analysis.stat import stat
from metaquantome.util.testutils import testfile, TTEST_SINFO
from metaquantome.util.constants import GO_TEST_DIR, EC_TEST_DIR


class TestFunctionalAnalysisExpand(unittest.TestCase):

    db = godb.GeneOntologyDb(GO_TEST_DIR, slim_down=True, overwrite=False)

    def testSingleInt(self):
        func=testfile('simple_func.tab')
        int=testfile('simple_int.tab')
        go_df = expand('f', sinfo='{"s1": ["int"]}', int_file=int, pep_colname_int='peptide',
                       pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=GO_TEST_DIR, func_file=func,
                       func_colname='go', ontology='go')
        self.assertEqual(go_df.loc["GO:0022610"]['int'], np.log2(200))
        self.assertEqual(go_df.loc["GO:0008152"]['int'], np.log2(100))

    def testMultipleInt(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        go_df = expand('f', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=int, pep_colname_int='peptide',
                       pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=GO_TEST_DIR, func_file=func,
                       func_colname='go', ontology='go')
        self.assertEqual(go_df.loc['GO:0008152']['int1'], np.log2(10))
        self.assertEqual(go_df.loc['GO:0022610']['int2'], np.log2(30))
        # missing values (zeros, nans, NA's, etc) are turned into NaN's
        self.assertTrue(np.isnan(go_df.loc['GO:0000003']['int3']))
        return go_df

    def testNopep(self):
        nopep=testfile('nopep.tab')
        go_df = expand('f', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=None, pep_colname_int='peptide',
                       pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=GO_TEST_DIR, func_colname='go',
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
        go_df = expand('f', sinfo=sinfo, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
                       pep_colname_tax='peptide', data_dir=GO_TEST_DIR, func_file=func, func_colname='go',
                       ontology='go', slim_down=True)
        # test that all go terms are in slim
        # load slim
        returned_gos = set(go_df['id'])
        # potential of unknown, so just drop that
        returned_gos.discard('unknown')
        self.assertTrue(returned_gos.issubset(self.db.goslim.keys()))

    def testCog(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        cog_df = expand('f', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=int, pep_colname_int='peptide',
                        pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=None, func_file=func,
                        func_colname='cog', ontology='cog')
        self.assertEqual(cog_df.loc["C"]['s1_mean'], np.log2((10+20+70)/3))
        self.assertEqual(cog_df.loc["N"]['int2'], np.log2(30))

    def testSimpleEc(self):
        func=testfile('simple_ec.tab')
        int=testfile('simple_int.tab')
        ec_df = expand('f', sinfo='{"s1": ["int"]}', int_file=int, pep_colname_int='peptide',
                       pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=EC_TEST_DIR, func_file=func,
                       func_colname='ec', ontology='ec')
        self.assertEqual(ec_df.loc["3.4.11.-"]['int'], np.log2(100))
        self.assertEqual(ec_df.loc["3.4.-.-"]['int'], np.log2(300))

    def testMultipleEc(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')
        ec_df = expand('f', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=int, pep_colname_int='peptide',
                       pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=EC_TEST_DIR, func_file=func,
                       func_colname='ec', ontology='ec')
        self.assertEqual(ec_df.loc['3.4.-.-']['int1'], np.log2(50))
        self.assertEqual(ec_df.loc['1.2.-.-']['int2'], np.log2(50))
        # missing values (zeros, nans, NA's, etc) are turned into NaN's
        self.assertTrue(np.isnan(ec_df.loc['1.2.-.-']['int3']))


class TestFunctionalAnalysisTest(unittest.TestCase):
    go_db = godb.GeneOntologyDb(GO_TEST_DIR, slim_down=True, overwrite=False)

    def testDA(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        expanded=testfile('go_expanded_ttest.tab')
        test_write=testfile('go_tested.tab')
        df_expd = expand('f', sinfo=TTEST_SINFO, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
                         pep_colname_tax='peptide', data_dir=GO_TEST_DIR, outfile=expanded, func_file=func,
                         func_colname='go', ontology='go')
        df_tst = stat(expanded, sinfo=TTEST_SINFO, paired=False, parametric=True, ontology='go', mode='f',
                      outfile=test_write)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(df_tst['p']['GO:0008152'] > 0.05)
        self.assertTrue(df_tst['p'][['GO:0022610','GO:0000003','GO:0032505']].le(0.05).all())

    def testCogTTest(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        expandfile=testfile('cog_ttest.tab')
        cog_df = expand('f', sinfo=TTEST_SINFO, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
                        pep_colname_tax='peptide', data_dir=None, outfile=expandfile, func_file=func,
                        func_colname='cog', ontology='cog')
        cog_tst = stat(expandfile, sinfo=TTEST_SINFO, paired=False, parametric=True, ontology='cog', mode='f',
                       outfile=None)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(cog_tst['p']['C'] > 0.05)
        self.assertTrue(cog_tst['p'][['N', 'D']].le(0.05).all())

    def testDiffAbundEc(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        expandfile=testfile('ec_ttest.tab')
        tested_file=testfile('ec_ttest_tested.tab')
        expand('f', sinfo=TTEST_SINFO, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
               pep_colname_tax='peptide', data_dir=EC_TEST_DIR, outfile=expandfile, func_file=func, func_colname='ec',
               ontology='ec')
        ec_tst = stat(expandfile, sinfo=TTEST_SINFO, paired=False, parametric=True, ontology='ec', mode='f',
                      outfile=tested_file)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(ec_tst['p']['3.4.11.-'] > 0.05)
        self.assertTrue(ec_tst['p'][['3.4.21.70', '1.2.-.-']].le(0.05).all())


if __name__ == '__main__':
    unittest.main()
