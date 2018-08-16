import unittest
from metaquant import taxonomy_database
from metaquant.definitions import DATA_DIR
import os
import metaquant.taxonomy_database as td
import shutil
from tests.testutils import testfile
import pandas as pd
from metaquant.definitions import DATA_DIR


class TestTaxonomyDatabase(unittest.TestCase):
    def testDownloadTaxonomy(self):
        # make tmp dir for testing download
        tmp_dir = os.path.join(DATA_DIR, 'tmp_test_tax_dwnld')
        os.mkdir(tmp_dir)
        ncbi = td.ncbi_database_handler(tmp_dir)
        lineage = ncbi.get_lineage(1919)
        shutil.rmtree(tmp_dir)
        self.assertTrue(1760 in lineage)

    def testMapToRanks(self):
        ncbi = taxonomy_database.ncbi_database_handler(DATA_DIR)
        # human genus
        # we don't want to get species back, because that's lower in the hierarchy
        human_ancestors = td.map_id_to_desired_ranks({'order', 'species', 'genus'}, 9605, ncbi)
        self.assertDictEqual({'order': 9443, 'genus': 9605}, human_ancestors)

    def testRealRankNamesMap(self):
        # take all taxids from unipept sample 7
        # just make sure that there aren't any ranks returned by unipept that are
        # not contained in FULL_TAXONOMY_TREE
        tax = testfile('unipept_sample7_taxonomy.tab')
        df = pd.read_csv(tax, sep='\t')
        ncbi = td.ncbi_database_handler(DATA_DIR)
        num_ids = td.convert_name_to_taxid(df['lca'], ncbi)
        for i in num_ids:
            td.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, i, ncbi)

    def testTranslator(self):
        ncbi = taxonomy_database.ncbi_database_handler(DATA_DIR)
        # has to be a list or series
        ecoli = [562]

        # returns list
        name = taxonomy_database.convert_taxid_to_name(ecoli, ncbi)
        self.assertEqual(name, ['Escherichia coli'])

    def testTranslateName2Id(self):
        ncbi = taxonomy_database.ncbi_database_handler(DATA_DIR)

        name= ['Brassicaceae']

        # returns list
        id = taxonomy_database.convert_name_to_taxid(name, ncbi)
        self.assertEqual(id, [3700])

    def testUnknownNames(self):
        # assigns the ncbi taxid for unassigned (32644) if a name is unknown

        ncbi = taxonomy_database.ncbi_database_handler(DATA_DIR)

        name = ['Random nonsense']

        id = taxonomy_database.convert_name_to_taxid(name, ncbi)
        self.assertEqual(id, [12908])

    def testRoot(self):
        ncbi = taxonomy_database.ncbi_database_handler(DATA_DIR)

        name = ['root']

        id = taxonomy_database.convert_name_to_taxid(name, ncbi)
        self.assertEqual(id, [1])

    def testUnipeptResults(self):
        # test a whole bunch of unipept name results
        ncbi = taxonomy_database.ncbi_database_handler(DATA_DIR)

        unipept_thaliana = os.path.join(DATA_DIR, 'test', 'unipept_a_thaliana_result.csv')
        with open(unipept_thaliana, mode = 'r') as f:
            f.readline() # ditch column header
            names = [elem.strip('\n') for elem in f.readlines()]
            ids = taxonomy_database.convert_name_to_taxid(names, ncbi)

        # assure no unknown names (i.e., unipept names are reasonably well supported)
        self.assertEqual(sum([id == 32644 for id in ids]), 0)



if __name__=='__main__':
    unittest.main()
