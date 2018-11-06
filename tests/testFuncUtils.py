from metaquant.databases import GeneOntologyDb as godb
from metaquant.util.utils import DATA_DIR
from metaquant.util import funcutils as fu
from metaquant.databases import EnzymeDb as ecdb

import unittest
import os
import pandas as pd


class TestFuncUtils(unittest.TestCase):
    GO_TEST_DIR = os.path.join(DATA_DIR, 'test', 'go_cache')  # downloaded 8/27/18
    go_db = godb.GeneOntologyDb(GO_TEST_DIR, slim_down=True)
    EC_TEST_DIR = os.path.join(DATA_DIR, 'test', 'ec_cache')  # downloaded 8/28/18
    ec_db = ecdb.EnzymeDb(EC_TEST_DIR)

    def testSplitFuncListGO(self):
        golist = 'go1,go2,go3'
        split = fu.split_func_list(golist, ',')
        exp_set = {'go1', 'go2', 'go3'}
        self.assertIsInstance(split, set)
        self.assertSetEqual(split, exp_set)

    def testMakeListNonRedundantGO(self):
        gochar = 'GO:0006915,GO:0008150'
        exp_result = 'GO:0006915'
        obs_result = fu.reduce_func(self.go_db, gochar, ',')
        self.assertEqual(exp_result, obs_result)

    def testUnipeptMakeListNonRedundantGO(self):
        gochar = 'GO:1903494,GO:1903496,GO:0042742,GO:0005794,GO:0005796,GO:0005576,GO:0042803,GO:0035375'
        obs_result = fu.reduce_func(self.go_db, gochar, ',')
        self.assertSetEqual(set(obs_result.split(',')), set(gochar.split(',')))

    def testMakeListNonRedundantEC(self):
        ecchar = '1.1.4.-,1.1.4.2'
        exp_result = '1.1.4.2'
        obs_result = fu.reduce_func(self.ec_db, ecchar, ',')
        self.assertEqual(exp_result, obs_result)

    def testMakeDfNonRedundantGO(self):
        godf = pd.DataFrame({'go': ['GO:0006915,GO:0008150', 'GO:0008150'],
                             'peptide': ['A', 'B'],
                             'intensity': [100, 200]}).sort_index(axis='columns')
        obs_result = fu.reduce_func_df(self.go_db, godf, 'go', ',')
        expected = pd.DataFrame({'go': ['GO:0006915', 'GO:0008150'],
                                 'peptide': ['A', 'B'],
                                 'intensity': [100, 200]}).sort_index(axis='columns')
        self.assertTrue(expected.equals(obs_result))

    def testMakeDfNonRedundantEC(self):
        ecdf = pd.DataFrame({'ec': ['1.1.4.-,1.1.4.2', '1.1.4.-,1.1.4.1'],
                             'peptide': ['A', 'B'],
                             'intensity': [100, 200]}).sort_index(axis='columns')
        obs_result = fu.reduce_func_df(self.ec_db, ecdf, 'ec', ',')
        expected = pd.DataFrame({'ec': ['1.1.4.2', '1.1.4.1'],
                                 'peptide': ['A', 'B'],
                                 'intensity': [100, 200]}).sort_index(axis='columns')
        self.assertTrue(expected.equals(obs_result))
