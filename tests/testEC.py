import unittest
from metaquant.runner import metaquant
from metaquant import ec
from tests.testutils import testfile
from metaquant.definitions import DATA_DIR
import os
import shutil
import numpy as np


class TestEC(unittest.TestCase):
    def testDownloadAndConversion(self):
        # make tmp dir for testing download
        tmp_dir = os.path.join(DATA_DIR, 'tmp_test_data_dwnld')
        os.mkdir(tmp_dir)
        enzyme_db = ec.enzyme_database_handler(tmp_dir, False)
        expected_contents = [os.path.join(tmp_dir, file)
                             for file in ['enzclass.txt', 'enzyme.dat', 'ec_id.json', 'enzclass.json']]
        for content in expected_contents:
            self.assertTrue(os.path.exists(content))
        shutil.rmtree(tmp_dir)

        # make sure parsed correctly
        # this is from enzyme.dat
        self.assertEqual(enzyme_db['1.2.3.4'], 'Oxalate oxidase.')

        # from enzclass.txt
        self.assertEqual(enzyme_db['6.1.-.-'], 'Forming carbon-oxygen bonds.')

    def testExpandEC(self):
        test_ec = '1.2.-.-'
        df = ec.expand_ec(test_ec)
        # should have two levels (the first two)
        self.assertEqual(df.loc[ec.LEVEL_NAMES[0]], '1.-.-.-')
        self.assertEqual(df.loc[ec.LEVEL_NAMES[1]], '1.2.-.-')
        self.assertTrue(ec.LEVEL_NAMES[2] not in df.index)

    def testSingleInt(self):
        func=testfile('simple_ec.tab')
        int=testfile('simple_int.tab')

        ec_df = metaquant('fn', sample_names={'s1': ['int']}, int_file=int, pep_colname='peptide', func_file=func,
                          ontology='ec', test=False, overwrite=False)
        # leaf of tree
        self.assertEqual(ec_df.loc['3.4.21.70']['int'], np.log2(200))

        # internal node - check that we are adding up the hierarchy
        self.assertEqual(ec_df.loc['3.4.-.-']['int'], np.log2(100+200))

    def testUnknownEC(self):
        func=testfile('unk_ec.tab')
        int=testfile('simple_int.tab')
        ec_df = metaquant('fn', sample_names={'s1': ['int']}, int_file=int, pep_colname='peptide', func_file=func,
                          ontology='ec', test=False, overwrite=False)
        self.assertEqual(ec_df.loc['1.50.10000.-']['description'], 'unknown_ec')



