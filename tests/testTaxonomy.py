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
                                   sample_names = {'samp1': ['intensity']})
        self.assertEqual(trey.query("id == 'Helicobacter pylori'")['samp1_mean'].values, 0.2)

    def testWrite(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_simple.tab')
        outfile = os.path.join(DATA_DIR, 'test', 'taxonomy_write_simple.tab')

        metaquant.metaquant(mode='tax', file=datafile,
                                   sample_names={'samp1': ['intensity']},
                                   outfile=outfile)
        written = pd.read_table(outfile)
        self.assertEqual(written.query("id == 'Clostridioides'")['samp1_mean'].values[0], 0.7)

    def testMultCols(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_test_multiple.tab')
        outfile = os.path.join(DATA_DIR, 'test', 'taxonomy_write_multiple.tab')

        samp_names = {'samp1': ['int1', 'int2'],
                      'samp2': ['int3', 'int4']}
        # without testing, make sure the different ranks sum to 1 in each group
        trey = metaquant.metaquant('tax',file=datafile,
                                   sample_names=samp_names,
                                   test=False,
                                   threshold=0)
        self.assertTrue(np.isclose(trey[['int1', 'int2', 'int3', 'int4', 'rank']].groupby(by="rank").sum(axis=0),1.0).all())

        # analysing, writing to file
        trey = metaquant.metaquant('tax',file=datafile,
                                   sample_names=samp_names,
                                   test=True,
                                   threshold=2,
                                   outfile=outfile)
        self.assertEqual(trey.query("rank == 'phylum' and id == 'Proteobacteria'")['int3'].values[0], 9000/10100)
        print(trey[trey.id == "Helicobacter pylori"]['log2fc_samp1_over_samp2'].values)

        # fold change
        expected = [np.log2(((2/10 + 1/2)/2)/((9000/10100 + 2/3)/2))]
        self.assertTrue(np.isclose(trey[trey.id == "Helicobacter pylori"]['log2fc_samp1_over_samp2'], expected))


if __name__ == '__main__':
    unittest.main()