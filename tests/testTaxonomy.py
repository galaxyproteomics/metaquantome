import unittest

import numpy as np
import pandas as pd

from src import taxonomy as tax


class TestTaxonomy(unittest.TestCase):
    def testSingleBasic(self):
        trey = tax.taxonomy_analysis('data/test/taxonomy_simple.tab', ['intensity'])
        self.assertEqual(trey.loc['pylori', 'intensity'], 0.2)

    def testFullDf(self):
        trey = tax.taxonomy_analysis('data/test/taxonomy_unipept_results.tabular', ['intensity'])
        eik = trey.query("rank == 'genus' and member == 'Eikenella'")['intensity'].values[0]
        self.assertAlmostEqual(eik, 0.0000337, places = 3)

    def testWrite(self):
        trey = tax.taxonomy_analysis('data/test/taxonomy_simple.tab', ['intensity'],
                                     outfile='data/test/taxonomy_write_simple.tab')
        written = pd.read_table('data/test/taxonomy_write_simple.tab')
        self.assertEqual(written.query("member == 'clostridium'")['intensity'].values[0], 0.7)

    def testMultCols(self):
        trey = tax.taxonomy_analysis('data/test/taxonomy_test_multiple.tab', ['int1', 'int2'], ['int3', 'int4'])
        self.assertEqual(trey.query("rank == 'phylum' and member == 'proteobacteria'")['int3'].values[0], 0.9)

    def testT(self):
        trey = tax.taxonomy_analysis('data/test/taxonomy_test_multiple.tab',
                                     ['int1','int2'],
                                     ['int3', 'int4'],
                                     test=True)
        expected = np.log2(((0.9 + 0.75)/2)/((0.2 + 0.5)/2))
        self.assertEqual(expected, trey[trey.member == "pylori"]['log2ratio_grp2_over_grp1'][0])

    def testReal(self):
        tax.taxonomy_analysis('data/test/taxonomy_unipept_3samps.tabular',
                              sample1_colnames=['int737NS', 'int852NS', 'int867NS'],
                              sample2_colnames=['int737WS', 'int852WS', 'int867WS'])


if __name__ == '__main__':
    unittest.main()