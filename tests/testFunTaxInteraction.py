import unittest
import metaquant
from definitions import DATA_DIR
import os
from src import stats
from tests.testutils import testfile
import numpy as np


class TestFunctionTaxInteraction(unittest.TestCase):
    def testReadAndDE(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        tax=testfile('multiple_tax.tab')

        ft_df = metaquant.metaquant('taxfn', sample_names={'s1': ['int1', 'int2', 'int3'],
                                                           's2': ['int4', 'int5', 'int6']}, int_file=int,
                                    func_file=func, tax_file=tax, ontology='cog', tax_colname='lca', test=True)

        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(ft_df.loc[(ft_df['taxon_name'] == 'Helicobacter pylori') & (ft_df['cog'] == 'C')]['corrected_p'].ge(0.05).all())
        self.assertTrue(ft_df.loc[ft_df['taxon_name'] != 'Helicobacter pylori']['corrected_p'].le(0.05).all())

    def testSumming(self):
        func=testfile('mult_ft_func.tab')
        int=testfile('int_ttest.tab')
        tax=testfile('mult_ft_tax.tab')

        ft_df = metaquant.metaquant('taxfn', sample_names={'s1': ['int1', 'int2', 'int3'],
                                                           's2': ['int4', 'int5', 'int6']}, int_file=int,
                                    func_file=func, tax_file=tax, ontology='cog', tax_colname='lca', test=True)

        # we expect the C-Helicobacter pylori pair to be the log of the mean intensity
        # {'int1': [10, 20, 1000],
        #  'int2': [20, 30, 1200],
        #  'int3': [15, 20, 900],
        #  'int4': [30, 3500, 12],
        #  'int5': [21, 2000, 13],
        #  'int6': [30, 3000, 10]}

        helico = ft_df.loc[(ft_df['taxon_name'] == 'Helicobacter pylori') & (ft_df['cog'] == 'C')]
        # test mean
        expected_s1_mean = np.log2((10 + 20 + 15 + 1000 + 1200 + 900)/3)
        self.assertTrue(helico['s1_mean'].eq(expected_s1_mean).all())

        # test log2fold change - is log2(s1_mean) - log2(s1_mean)
        expected_log1o2 = np.log2(np.mean([10, 1000, 20, 1200, 15, 900])) - np.log2(np.mean([30, 21, 30, 12, 13, 10]))
        self.assertTrue(ft_df.loc[(ft_df['taxon_name'] == 'Helicobacter pylori') & (ft_df['cog'] == 'C')]['log2fc_s1_over_s2'].eq(expected_log1o2).all())




if __name__=='__main__':
    unittest.main()