import unittest
import pandas as pd
import numpy as np

import metaquantome.util.utils
from metaquantome.databases import GeneOntologyDb as godb
from metaquantome.util import utils as utils
from metaquantome.databases import EnzymeDb as ecdb
from metaquantome.databases import NCBITaxonomyDb as ncbi
from metaquantome.util.utils import TEST_DIR


class TestUtils(unittest.TestCase):
    def testSafeCastToList(self):
        string = "hello"
        self.assertEqual(utils.safe_cast_to_list(string), [string])
        self.assertEqual(utils.safe_cast_to_list([string]), [string])

    def testSniffTaxNames(self):
        num_df = pd.DataFrame({'tax': ['122', '133', '999', '444']})
        self.assertFalse(utils.sniff_tax_names(num_df, 'tax'))

        nam_df = pd.DataFrame({'tax': ['imataxon', 'imanothertaxon', 'third']})
        self.assertTrue(utils.sniff_tax_names(nam_df, 'tax'))

    def testFilter(self):
        ncbi_db = ncbi.NCBITaxonomyDb(TEST_DIR)
        tax_df = pd.DataFrame({'tax': [np.nan, 210, '9999999']},
                              index=[0, 1, 2])
        filt_df = utils.filter_df(ncbi_db, 'tax', tax_df)
        self.assertIn(1, filt_df.index)
        self.assertNotIn(0, filt_df.index)
        self.assertNotIn(2, filt_df.index)


class TestReduceFuncDf(unittest.TestCase):
    go_db = godb.GeneOntologyDb(TEST_DIR, slim_down=True)
    ec_db = ecdb.EnzymeDb(TEST_DIR)

    def testMakeListNonRedundantGO(self):
        gochar = 'GO:0006915,GO:0008150'
        exp_result = 'GO:0006915'
        obs_result = metaquantome.util.utils.reduce_func(self.go_db, gochar, ',')
        self.assertEqual(exp_result, obs_result)

    def testUnipeptMakeListNonRedundantGO(self):
        gochar = 'GO:1903494,GO:1903496,GO:0042742,GO:0005794,GO:0005796,GO:0005576,GO:0042803,GO:0035375'
        obs_result = metaquantome.util.utils.reduce_func(self.go_db, gochar, ',')
        self.assertSetEqual(set(obs_result.split(',')), set(gochar.split(',')))

    def testMakeListNonRedundantEC(self):
        ecchar = '1.1.4.-,1.1.4.2'
        exp_result = '1.1.4.2'
        obs_result = metaquantome.util.utils.reduce_func(self.ec_db, ecchar, ',')
        self.assertEqual(exp_result, obs_result)

    def testMakeDfNonRedundantGO(self):
        godf = pd.DataFrame({'go': ['GO:0006915,GO:0008150', 'GO:0008150'],
                             'peptide': ['A', 'B'],
                             'intensity': [100, 200]}).sort_index(axis='columns')
        obs_result = metaquantome.util.utils.reduce_func_df(self.go_db, godf, 'go', ',')
        expected = pd.DataFrame({'go': ['GO:0006915', 'GO:0008150'],
                                 'peptide': ['A', 'B'],
                                 'intensity': [100, 200]}).sort_index(axis='columns')
        self.assertTrue(expected.equals(obs_result))

    def testMakeDfNonRedundantEC(self):
        ecdf = pd.DataFrame({'ec': ['1.1.4.-,1.1.4.2', '1.1.4.-,1.1.4.1'],
                             'peptide': ['A', 'B'],
                             'intensity': [100, 200]}).sort_index(axis='columns')
        obs_result = metaquantome.util.utils.reduce_func_df(self.ec_db, ecdf, 'ec', ',')
        expected = pd.DataFrame({'ec': ['1.1.4.2', '1.1.4.1'],
                                 'peptide': ['A', 'B'],
                                 'intensity': [100, 200]}).sort_index(axis='columns')
        self.assertTrue(expected.equals(obs_result))

if __name__ == "__main__":
    unittest.main()