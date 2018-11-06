import unittest
from metaquant.runner import metaquant_runner
from metaquant.analysis.expand import expand
from tests.testutils import testfile
import tests.testutils as tu
import numpy as np


class TestFunctionTaxInteraction(unittest.TestCase):
    def testReadAndDE(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')
        tax=testfile('multiple_tax.tab')
        ft_df = expand(mode='taxfn', samps=tu.TTEST_SINFO, int_file=int, func_colname='go', func_file=func,
                       tax_file=tax, ontology='go', tax_colname='lca',
                       slim_down=True)

        print(ft_df)
        # make sure calculated mean is accurate
        self.assertTrue(ft_df.loc[(ft_df['taxon_name'] == 'Helicobacter pylori') &
                                  (ft_df['go'] == 'GO:0008152')]['s1_mean'].ge(0.05).all())



if __name__=='__main__':
    unittest.main()