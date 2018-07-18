import unittest
from metaquant.runner import metaquant
from metaquant import ec
from tests.testutils import testfile
from metaquant.definitions import DATA_DIR
import os
import shutil
import numpy as np
import pandas as pd


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
                          ontology='ec', test=False, overwrite=False, func_colname='ec')
        # leaf of tree
        self.assertEqual(ec_df.loc['3.4.21.70']['int'], np.log2(200))

        # internal node - check that we are adding up the hierarchy
        self.assertEqual(ec_df.loc['3.4.-.-']['int'], np.log2(100+200))

    def testUnknownEC(self):
        func=testfile('unk_ec.tab')
        int=testfile('simple_int.tab')
        ec_df = metaquant('fn', sample_names={'s1': ['int']}, int_file=int, pep_colname='peptide', func_file=func,
                          ontology='ec', test=False, overwrite=False, func_colname='ec')
        self.assertEqual(ec_df.loc['1.50.10000.-']['description'], 'unknown_ec')

    def testExpandList(self):
        peptide = ['AAYEEAEHAAK', 'AGVTK', 'FAKE']
        df = pd.DataFrame({'EC': ['1.11.1.1,1.-.-.-,1.14.13.81', '2.3.1.234,1.2.1.-,4.2.1.9', '']},
                          index=peptide)
        print(ec.split_ec_list(df, 'EC'))

    def testRealEC(self):
        func=testfile('unipept_sample7_functional_clean.tab')
        int=testfile('unipept_sample7_int_clean.tab')
        ec_df = metaquant('fn', sample_names={'s1': ['int']}, int_file=int, pep_colname='peptide', func_file=func,
                          ontology='ec', test=False, overwrite=False, func_colname='EC')
        # make sure that all of the 1s have been filtered out
        self.assertEqual(ec_df.query('id == "1.-.-.-"').size, 0)


