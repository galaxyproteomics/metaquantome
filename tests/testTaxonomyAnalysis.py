import unittest
import numpy as np
import pandas as pd
from metaquant.runner import metaquant_runner
from tests.testutils import testfile, TTEST_SINFO


class TestTaxonomyAnalysis(unittest.TestCase):

    def testSingleBasic(self):
        tax = testfile('simple_tax.tab')
        int = testfile('simple_int.tab')
        tax_df = metaquant_runner('tax', sinfo='{"s1": ["int"]}', int_file=int, pep_colname='peptide', tax_file=tax,
                                  tax_colname='lca', test=False)
        self.assertEqual(tax_df.query("taxon_name == 'Helicobacter pylori'")['int'].values, np.log2(100))

    def testWrite(self):
        tax = testfile('simple_tax.tab')
        int = testfile('simple_int.tab')
        out = testfile('taxonomy_write_simple.tab')
        metaquant_runner(mode='tax', sinfo='{"samp1": ["int"]}', int_file=int, tax_file=tax, tax_colname='lca', outfile=out)
        written = pd.read_table(out)
        self.assertAlmostEqual(written.query("taxon_name == 'Clostridioides difficile'")['samp1_mean'].values[0], np.log2(200))

    def testMultCols(self):
        tax=testfile('multiple_tax.tab')
        int=testfile('multiple_int.tab')
        tax_df = metaquant_runner('tax', sinfo='{"s1": ["int1", "int2", "int3"]}', int_file=int, tax_file=tax, tax_colname='lca')
        self.assertEqual(tax_df.query("rank == 'phylum' and taxon_name == 'Proteobacteria'")['int3'].values[0], np.log2(70))

    def testTaxTTests(self):
        tax=testfile('multiple_tax.tab')
        int=testfile('int_ttest.tab')
        tax_df = metaquant_runner('tax', sinfo=TTEST_SINFO, int_file=int, tax_file=tax, tax_colname='lca', test=True,
                                  paired=False, parametric=False)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(tax_df['p'][210] > 0.05)
        self.assertTrue(tax_df['p'][[1496,1870884]].le(0.05).all())
        # also, make sure firmicutes phylum is sum of c difficile and clostridiaceae, divided by all phyla
        self.assertEqual(tax_df['int1'][1239], np.log2(1020))


if __name__ == '__main__':
    unittest.main()