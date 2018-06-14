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

        ft_df = metaquant.metaquant('taxfn', func_file=func,
                                    int_file=int,
                                    tax_file=tax,
                                    tax_colname='lca',
                                    ontology='cog',
                                    sample_names={'s1': ['int1', 'int2', 'int3'],
                                                  's2': ['int4', 'int5', 'int6']},
                                    test=True)

        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(ft_df['corrected_p']['C-Helicobacter pylori'] > 0.05)
        self.assertTrue(ft_df['corrected_p'][['D-Clostridioides','N-Clostridioides difficile']].le(0.05).all())

    def testSumming(self):
        func=testfile('mult_ft_func.tab')
        int=testfile('int_ttest.tab')
        tax=testfile('mult_ft_tax.tab')

        ft_df = metaquant.metaquant('taxfn', func_file=func,
                                    int_file=int,
                                    tax_file=tax,
                                    tax_colname='lca',
                                    ontology='cog',
                                    sample_names={'s1': ['int1', 'int2', 'int3'],
                                                  's2': ['int4', 'int5', 'int6']},
                                    test=True)

        # we expect the C-Helicobacter pylori pair to be the log of the mean intensity
        # {'int1': [10, 20, 1000],
        #  'int2': [20, 30, 1200],
        #  'int3': [15, 20, 900],
        #  'int4': [30, 3500, 12],
        #  'int5': [21, 2000, 13],
        #  'int6': [30, 3000, 10]}

        # test mean
        expected_s1_mean = np.log2((10 + 20 + 15 + 1000 + 1200 + 900)/3)
        self.assertEqual(ft_df['s1_mean']['C-Helicobacter pylori'], expected_s1_mean)

        # test log2fold change - is log2(s1_mean) - log2(s1_mean)
        expected_log1o2 = np.log2(np.mean([10, 1000, 20, 1200, 15, 900])) - np.log2(np.mean([30, 21, 30, 12, 13, 10]))
        self.assertEqual(ft_df['log2fc_s1_over_s2']['C-Helicobacter pylori'], expected_log1o2)




if __name__=='__main__':
    unittest.main()