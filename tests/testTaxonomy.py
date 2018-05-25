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
        self.assertEqual(trey.loc['pylori', 'samp1_mean'], 0.2)

    def testFullDf(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_unipept_results.tabular')
        trey = metaquant.metaquant('tax', file=datafile, sample_names={'samp1' : ['intensity']})
        # the column named intensity is kept by default
        eik = trey.query("rank == 'genus' and id == 'Eikenella'")['intensity'].values[0]
        self.assertAlmostEqual(eik, 0.0000337, places=3)

    def testWrite(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_simple.tab')
        outfile = os.path.join(DATA_DIR, 'test', 'taxonomy_write_simple.tab')

        metaquant.metaquant(mode='tax', file=datafile,
                                   sample_names={'samp1': ['intensity']},
                                   outfile=outfile)
        written = pd.read_table(outfile)
        self.assertEqual(written.query("id == 'clostridium'")['samp1_mean'].values[0], 0.7)

    def testMultCols(self):
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_test_multiple.tab')
        outfile = os.path.join(DATA_DIR, 'test', 'taxonomy_write_multiple.tab')

        samp_names = {'samp1': ['int1', 'int2'],
                      'samp2': ['int3', 'int4']}
        # without testing, make sure the different ranks sum to 1 in each group
        trey = metaquant.metaquant('tax',file=datafile,
                                   sample_names=samp_names,
                                   test=False,
                                   threshold=2)
        self.assertTrue(np.isclose(trey[['int1', 'int2', 'int3', 'int4', 'rank']].groupby(by="rank").sum(axis=0),1.0).all())

        # analysing, writing to file
        trey = metaquant.metaquant('tax',file=datafile,
                                   sample_names=samp_names,
                                   test=True,
                                   threshold=2,
                                   outfile=outfile)
        self.assertEqual(trey.query("rank == 'phylum' and id == 'proteobacteria'")['int3'].values[0], 9/11)

        # fold change
        expected = np.log2(((2/10 + 1/2)/2)/((9/11 + 2/3)/2))
        self.assertEqual(expected, trey[trey.id == "pylori"]['log2fc_samp1_over_samp2'][0])

    def testTaxFiltering(self):
        """
        test that only tax ranks that are known at least once are kept in dataset
        """
        datafile = os.path.join(DATA_DIR, 'test', 'taxonomy_unipept_small_3samps.tabular')
        sample_names={'NS': ['int737NS', 'int852NS', 'int867NS'],
                      'WS': ['int737WS', 'int852WS', 'int867WS']}
        samps = metaquant.common.SampleGroups(sample_names)

        df = metaquant.common.read_data_table(datafile, samps, "peptide")
        expected_cols = {'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}
        true_cols = set(df.drop(samps.all_intcols, axis = 1))
        self.assertEqual(expected_cols, true_cols)


if __name__ == '__main__':
    unittest.main()