import unittest
from metaquant.runner import metaquant_runner
from tests.testutils import testfile
import tests.testutils as tu
import numpy as np


class TestFunctionTaxInteraction(unittest.TestCase):
    def testReadAndDE(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        tax=testfile('multiple_tax.tab')
        ft_df = metaquant_runner(mode='taxfn', sinfo=tu.TTEST_SINFO, int_file=int, func_colname='cog', func_file=func,
                                 tax_file=tax, ontology='cog', tax_colname='lca', test=True, parametric=True)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(ft_df.loc[(ft_df['taxon_name'] == 'Helicobacter pylori') & (ft_df['cog'] == 'C')]['corrected_p'].ge(0.05).all())
        self.assertTrue(ft_df.loc[ft_df['taxon_name'] != 'Helicobacter pylori']['corrected_p'].le(0.05).all())

    def testSumming(self):
        func=testfile('mult_ft_func.tab')
        int=testfile('int_ttest.tab')
        tax=testfile('mult_ft_tax.tab')
        ft_df = metaquant_runner('taxfn', sinfo=tu.TTEST_SINFO, int_file=int, func_colname='cog', func_file=func, tax_file=tax,
                                 ontology='cog', tax_colname='lca', test=True, parametric=True)
        # we expect the C-Helicobacter pylori pair to be the log of the mean intensity
        # {'int1': [12, 20, 1000],
        #  'int2': [20, 30, 1200],
        #  'int3': [15, 20, 900],
        #  'int4': [12, 3500, 12],
        #  'int5': [21, 2000, 13],
        #  'int6': [10, 3000, 10]}
        helico = ft_df.loc[(ft_df['taxon_name'] == 'Helicobacter pylori') & (ft_df['cog'] == 'C')]

        # test mean
        expected_s1_mean = np.log2((12 + 20 + 15 + 1000 + 1200 + 900)/3)
        self.assertTrue(helico['s1_mean'].eq(expected_s1_mean).all())

        # test log2fold change - is log2(s1_mean) - log2(s1_mean)
        expected_log1o2 = np.log2(np.mean([12, 1000, 20, 1200, 15, 900])) - np.log2(np.mean([12, 21, 10, 12, 13, 10]))
        self.assertTrue(np.isclose(
            ft_df.loc[(ft_df['taxon_name'] == 'Helicobacter pylori') & (ft_df['cog'] == 'C')]['log2fc_s1_over_s2'],
            expected_log1o2)
        )


if __name__=='__main__':
    unittest.main()