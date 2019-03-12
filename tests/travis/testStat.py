import unittest
import numpy as np
import pandas as pd

from metaquantome.databases import GeneOntologyDb as godb
import metaquantome.modules.expand as expand
from metaquantome.classes.SampleGroups import SampleGroups
import metaquantome.modules.stat as stat
from metaquantome.util import testutils as tu
from metaquantome.util.testutils import testfile, TTEST_SINFO
from metaquantome.util.utils import TEST_DIR


class TestStatUtils(unittest.TestCase):
    def testFoldChange(self):
        df = pd.DataFrame({'s1_1': [4, 4],
                           's1_2': [2, 2],
                           's2_1': [5, 10],
                           's2_2': [7, 16]})
        samps = SampleGroups('{"s1": ["s1_1", "s1_2"], "s2": ["s2_1", "s2_2"]}')

        means = expand.calc_means(df, samps)
        fc = stat.log2_fold_change(means, samps)
        self.assertTrue(fc['log2fc_s1_over_s2'].equals(pd.Series({0: np.log2(3/6), 1: np.log2(3/13)})))


class TestFunctionTaxInteraction(unittest.TestCase):
    def testReadAndDE(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        '''intensity:
        peptide	int1	int2	int3	int4	int5	int6
        A	12	20	15	12	21	10
        B	20	30	20	3500	2000	3000
        C	1000	1200	900	12	13	10
        '''
        # todo - add test for non-slim

        tax=testfile('multiple_tax.tab')
        ft_out=testfile('ft_out.tab')
        ft_df = expand.expand(mode='ft', sinfo=tu.TTEST_SINFO, int_file=int, pep_colname_int='peptide',
                              pep_colname_func='peptide', pep_colname_tax='peptide', data_dir=TEST_DIR, outfile=ft_out,
                              func_file=func, func_colname='go', ontology='go', slim_down=True, tax_file=tax,
                              tax_colname='lca')
        # make sure calculated mean is accurate
        # b and c both map to 8150
        exp_s1_mean = np.log2(((20 + 1000) + (1200 + 30) + (900 + 20))/3)
        obtained_mean = ft_df.loc[(ft_df['taxon_name'] == 'Clostridioides') &
                                  (ft_df['go'] == 'GO:0008150'), 's1_mean'][0]

        self.assertEqual(exp_s1_mean, obtained_mean)


class TestFunctionalAnalysisTest(unittest.TestCase):
    go_db = godb.GeneOntologyDb(TEST_DIR, slim_down=True)

    def testDA(self):
        func = testfile('multiple_func.tab')
        int = testfile('int_ttest.tab')
        expanded = testfile('go_expanded_ttest.tab')
        test_write = testfile('go_tested.tab')
        df_expd = expand.expand('f', sinfo=TTEST_SINFO, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
                                pep_colname_tax='peptide', data_dir=TEST_DIR, outfile=expanded, func_file=func,
                                func_colname='go', ontology='go')
        df_tst = stat.stat(expanded, sinfo=TTEST_SINFO, paired=False, parametric=True, ontology='go', mode='f',
                           outfile=test_write)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(df_tst['p']['GO:0008152'] > 0.05)
        self.assertTrue(df_tst['p'][['GO:0022610','GO:0000003','GO:0032505']].le(0.05).all())

    def testCogTTest(self):
        func = testfile('multiple_func.tab')
        int = testfile('int_ttest.tab')
        expandfile = testfile('cog_ttest.tab')
        cog_df = expand.expand('f', sinfo=TTEST_SINFO, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
                               pep_colname_tax='peptide', outfile=expandfile, func_file=func, func_colname='cog',
                               ontology='cog')
        cog_tst = stat.stat(expandfile, sinfo=TTEST_SINFO, paired=False, parametric=True, ontology='cog', mode='f',
                            outfile=None)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(cog_tst['p']['C'] > 0.05)
        self.assertTrue(cog_tst['p'][['N', 'D']].le(0.05).all())

    def testDiffAbundEc(self):
        func = testfile('multiple_func.tab')
        int = testfile('int_ttest.tab')
        expandfile = testfile('ec_ttest.tab')
        tested_file = testfile('ec_ttest_tested.tab')
        expand.expand('f', sinfo=TTEST_SINFO, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
                      pep_colname_tax='peptide', data_dir=TEST_DIR, outfile=expandfile, func_file=func, func_colname='ec',
                      ontology='ec')
        ec_tst = stat.stat(expandfile, sinfo=TTEST_SINFO, paired=False, parametric=True, ontology='ec', mode='f',
                           outfile=tested_file)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(ec_tst['p']['3.4.11.-'] > 0.05)
        self.assertTrue(ec_tst['p'][['3.4.21.70', '1.2.-.-']].le(0.05).all())


class TestTaxonomyAnalysisTest(unittest.TestCase):
    def testTaxTTests(self):
        tax = testfile('multiple_tax.tab')
        int = testfile('int_ttest.tab')
        expanded = testfile('expand_taxttest.tab')
        tax_df = expand.expand('t', sinfo=TTEST_SINFO, int_file=int, pep_colname_int='peptide', pep_colname_func='peptide',
                               pep_colname_tax='peptide', data_dir=TEST_DIR, outfile=expanded, tax_file=tax, tax_colname='lca')
        tax_tst = stat.stat(expanded, sinfo=TTEST_SINFO, paired=False, parametric=False, ontology=None, mode=None,
                            outfile=None)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(tax_tst['p'][210] > 0.05)
        self.assertTrue(tax_tst['p'][[1496,1870884]].le(0.05).all())
        # also, make sure firmicutes phylum is sum of c difficile and clostridiaceae
        self.assertEqual(tax_tst['int1'][1239], np.log2(1020))


if __name__ == '__main__':
    unittest.main()
