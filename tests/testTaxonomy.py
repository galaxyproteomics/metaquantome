import unittest
import numpy as np
import pandas as pd
import metaquant
import os
from definitions import DATA_DIR


class TestTaxonomy(unittest.TestCase):
    def testSingleBasic(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_simple.tab')
        trey = metaquant.metaquant('tax', file=datafile,
                                   sample1_colnames = 'intensity')
        self.assertEqual(trey.loc['pylori', 'intensity'], 0.2)

    def testFullDf(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_unipept_results.tabular')
        trey = metaquant.metaquant('tax', file=datafile, sample1_colnames='intensity')
        eik = trey.query("rank == 'genus' and member == 'Eikenella'")['intensity'].values[0]
        self.assertAlmostEqual(eik, 0.0000337, places = 3)

    def testWrite(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_simple.tab')
        outfile = os.path.join(DATA_DIR, 'test', 'taxonomy_write_simple.tab')

        metaquant.metaquant(mode='tax', file=datafile,
                                   sample1_colnames='intensity',
                                   outfile=outfile)
        written = pd.read_table(outfile)
        self.assertEqual(written.query("member == 'clostridium'")['intensity'].values[0], 0.7)

    def testMultCols(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_test_multiple.tab')
        outfile = os.path.join(DATA_DIR, 'test', 'taxonomy_write_multiple.tab')

        # analysing, writing to file
        trey = metaquant.metaquant('tax',file=datafile,
                                   sample1_colnames=['int1', 'int2'],
                                   sample2_colnames=['int3', 'int4'],
                                   test=True,
                                   threshold=2)
        self.assertEqual(trey.query("rank == 'phylum' and member == 'proteobacteria'")['int3'].values[0], 0.9)

        # fold change - this is assuming that clostridia gets filtered out
        expected = np.log2(((0.9 + 0.75)/2)/((2/9 + 5/7)/2))
        self.assertEqual(expected, trey[trey.member == "pylori"]['log2ratio_2over1'][0])


if __name__ == '__main__':
    unittest.main()