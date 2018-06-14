import unittest
import numpy as np
import pandas as pd
import metaquant
import os
from definitions import DATA_DIR
from tests.testutils import testfile

class TestTaxonomy(unittest.TestCase):
    def testSingleBasic(self):
        tax = testfile('simple_tax.tab')
        int = testfile('simple_int.tab')
        tax_df = metaquant.metaquant('tax', pep_colname='peptide',
                                    tax_file=tax, int_file=int,
                                    tax_colname='lca',
                                    sample_names={'s1': ['int']}, test=False)
        self.assertEqual(tax_df.query("id == 'Helicobacter pylori'")['int'].values, np.log2(100))

    def testWrite(self):
        tax = testfile('simple_tax.tab')
        int = testfile('simple_int.tab')
        out = testfile('taxonomy_write_simple.tab')

        metaquant.metaquant(mode='tax', int_file=int, tax_file=tax,
                            tax_colname='lca',
                            sample_names={'samp1': ['int']},
                            outfile=out)

        written = pd.read_table(out)
        self.assertAlmostEqual(written.query("id == 'Clostridioides'")['samp1_mean'].values[0], np.log2(200))

    def testMultCols(self):
        tax=testfile('multiple_tax.tab')
        int=testfile('multiple_int.tab')

        tax_df = metaquant.metaquant('tax',
                                     tax_file=tax,
                                     tax_colname='lca',
                                     int_file=int,
                                     sample_names={'s1': ['int1', 'int2', 'int3']})

        self.assertEqual(tax_df.query("rank == 'phylum' and id == 'Proteobacteria'")['int3'].values[0], np.log2(70))

    def testTaxTTests(self):
        tax=testfile('multiple_tax.tab')
        int=testfile('int_ttest.tab')

        tax_df = metaquant.metaquant('tax',
                                     tax_file=tax,
                                     tax_colname='lca',
                                     int_file=int,
                                     sample_names={'s1': ['int1', 'int2', 'int3'],
                                                   's2': ['int4', 'int5', 'int6']},
                                     test=True)

        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(tax_df['corrected_p'][210] > 0.05)
        self.assertTrue(tax_df['corrected_p'][[1496,1870884]].le(0.05).all())

        # also, make sure firmicutes phylum is sum of c difficile and clostridiaceae, divided by all phyla
        self.assertEqual(tax_df['int1'][1239], np.log2(1020))


if __name__ == '__main__':
    unittest.main()