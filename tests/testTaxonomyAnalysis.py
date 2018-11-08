import unittest
import numpy as np
import pandas as pd
from tests.testutils import testfile, TTEST_SINFO

from metaquantome.analysis.expand import expand
from metaquantome.analysis.stat import stat


class TestTaxonomyAnalysisExpand(unittest.TestCase):

    def testSingleBasic(self):
        tax = testfile('simple_tax.tab')
        int = testfile('simple_int.tab')
        tax_df = expand('tax', samps='{"s1": ["int"]}', int_file=int, pep_colname='peptide', tax_file=tax,
                        tax_colname='lca')
        self.assertEqual(tax_df.query("taxon_name == 'Helicobacter pylori'")['int'].values, np.log2(100))

    def testWrite(self):
        tax = testfile('simple_tax.tab')
        int = testfile('simple_int.tab')
        out = testfile('taxonomy_write_simple.tab')
        df = expand(mode='tax', samps='{"samp1": ["int"]}', int_file=int, outfile=out, tax_file=tax, tax_colname='lca')
        written = pd.read_table(out)
        self.assertAlmostEqual(written.query("taxon_name == 'Clostridioides difficile'")['samp1_mean'].values[0], np.log2(200))

    def testMultCols(self):
        tax=testfile('multiple_tax.tab')
        int=testfile('multiple_int.tab')
        tax_df = expand('tax', samps='{"s1": ["int1", "int2", "int3"]}', int_file=int, tax_file=tax, tax_colname='lca')
        self.assertEqual(tax_df.query("rank == 'phylum' and taxon_name == 'Proteobacteria'")['int3'].values[0], np.log2(70))

    def testNopep(self):
        nopep=testfile('nopep.tab')
        tax_df = expand('tax', samps='{"s1": ["int1", "int2", "int3"]}', tax_colname='lca', nopep=True,
                        nopep_file=nopep)
        self.assertEqual(tax_df.query("rank == 'phylum' and taxon_name == 'Proteobacteria'")['int3'].values[0],
                         np.log2(70))

    def testParentIntensityHigher(self):
        """
        make sure that parents always have higher intensity than children
        """
        tax=testfile('test_root_sum_uni.tab')
        int=testfile('test_root_sum_int.tab')
        tax_df = expand('tax', samps='{"A": ["int"]}', int_file=int, pep_colname='peptide', tax_file=tax,
                        tax_colname='taxon_id')
        # filter to phylum and below
        tax_df_filt = tax_df[(tax_df["rank"] != 'no rank') & (tax_df["rank"] != 'superkingdom')]
        # firmicutes phylum should be highest
        ints = tax_df_filt['int']
        self.assertEqual(ints.max(), ints[1239])
        # strep genus intensity should be greater than or equal to that of strep species
        self.assertGreaterEqual(ints[1301], ints[1302])
        self.assertGreaterEqual(ints[1301], ints[1305])


class TestTaxonomyAnalysisTest(unittest.TestCase):
    def testTaxTTests(self):
        tax=testfile('multiple_tax.tab')
        int=testfile('int_ttest.tab')
        expanded=testfile('expand_taxttest.tab')
        tax_df = expand('tax', samps=TTEST_SINFO, int_file=int, tax_file=tax, tax_colname='lca',
                        outfile=expanded)
        tax_tst = stat(expanded, samps=TTEST_SINFO, paired=False, parametric=False, ontology=None, mode=None,
                       outfile=None)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(tax_tst['p'][210] > 0.05)
        self.assertTrue(tax_tst['p'][[1496,1870884]].le(0.05).all())
        # also, make sure firmicutes phylum is sum of c difficile and clostridiaceae
        self.assertEqual(tax_tst['int1'][1239], np.log2(1020))


if __name__ == '__main__':
    unittest.main()