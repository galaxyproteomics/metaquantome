import unittest
import numpy as np

from metaquantome.analysis.expand import expand
from metaquantome.util.testutils import testfile
import metaquantome.util.testutils as tu
from metaquantome.util.constants import GO_TEST_DIR, TAX_TEST_DIR


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
        ft_df = expand(mode='taxfn', sinfo=tu.TTEST_SINFO, int_file=int, func_file=func, func_colname='go',
                       ontology='go', slim_down=True, tax_file=tax, tax_colname='lca',
                       ft_func_data_dir=GO_TEST_DIR, ft_tax_data_dir=TAX_TEST_DIR,
                       outfile=ft_out)
        # make sure calculated mean is accurate
        # b and c both map to 8150
        exp_s1_mean = np.log2(((20 + 1000) + (1200 + 30) + (900 + 20))/3)
        obtained_mean = ft_df.loc[(ft_df['taxon_name'] == 'Clostridioides') &
                                  (ft_df['go'] == 'GO:0008150'), 's1_mean'][0]

        self.assertEqual(exp_s1_mean, obtained_mean)


if __name__=='__main__':
    unittest.main()