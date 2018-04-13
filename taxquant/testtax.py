import unittest
import taxquant.tax as tax
import pandas as pd


class TestToyExample(unittest.TestCase):
    def testRelAbund(self):
        trey = tax.TreeOfLife('test/test.tab', ['intensity'])
        self.assertEqual(list(trey.rel_abundance_rank('genus')['intensity']), [0.7, 0.2, 0.1])


class TestUnipeptResults(unittest.TestCase):
    def testUnipeptIn(self):
        trey = tax.TreeOfLife('test/unipept_results.tabular', ['intensity'])
        all = trey.rel_abundance_all_ranks().query("rank == 'superkingdom' and member == 'Eukaryota'")['intensity'].values[0]
        self.assertAlmostEqual(all, 0.755, places = 3)


class TestAllRelAbund(unittest.TestCase):
    def testFullDf(self):
        trey = tax.TreeOfLife('test/unipept_results_metaproteomics.tabular', ['intensity'])
        eik = trey.rel_abundance_rank('genus').query("member == 'Eikenella'")['intensity'].values[0]
        self.assertAlmostEqual(eik, 0.0000337, places = 3)


class TestWrite(unittest.TestCase):
    def testRelAbund(self):
        trey = tax.TreeOfLife('test/test.tab', ['intensity'])
        trey.write_out('test/out.tab')
        written = pd.read_table('test/out.tab')
        self.assertEqual(written.query("member == 'clostridium'")['intensity'].values[0], 0.7)


class TestMultCols(unittest.TestCase):
    def testMultCols(self):
        trey = tax.TreeOfLife('test/test_multiple.tab', ['int1','int2'], ['int3', 'int4'])
        abund = trey.rel_abundance_all_ranks()
        self.assertEqual(abund.query("rank == 'phylum' and member == 'proteobacteria'")['int3'].values[0], 0.9)

class TestTTest(unittest.TestCase):
    def testT(self):
        trey = tax.TreeOfLife('test/test_multiple.tab', ['int1','int2'], ['int3', 'int4'])
        trey.test_norm_intensity()


if __name__ == '__main__':
    unittest.main()