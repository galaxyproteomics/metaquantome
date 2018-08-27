import unittest
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
        try:
            ncbi = td.NCBITaxonomyDb(tmp_dir)
            lineage = ncbi.get_ancestors(1919)
            self.assertTrue(1760 in lineage)
        finally:
            shutil.rmtree(tmp_dir)

    def testMapIdToDesiredRanks(self):
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        # human genus
        # we don't want to get species back, because that's lower in the hierarchy
        human_ancestors = ncbi.map_id_to_desired_ranks({'order', 'species', 'genus'}, 9605)
        self.assertDictEqual({'order': 9443, 'genus': 9605}, human_ancestors)

    def testMapIdToDesiredRanks_Unipept(self):
        # take all taxids from unipept sample 7
        # just make sure that there aren't any ranks returned by unipept that are
        # not contained in FULL_TAXONOMY_TREE
        tax = testfile('unipept_sample7_taxonomy.tab')
        df = pd.read_csv(tax, sep='\t')
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        num_ids = ncbi.convert_name_to_taxid(df['lca'])
        for i in num_ids:
            ncbi.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, i)

    def testExpandSampleTaxonomy(self):
        # test that expanded sample set contains all desired ancestors
        ncbi = td.NCBITaxonomyDb(DATA_DIR)

        # expect these taxids and all ancestors within BASIC_TAXONOMY_TREE
        sample_set = {9605, 562, 1807140}
        ecoli = set(ncbi.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 562).values())
        homo = set(ncbi.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 9605).values())
        acidithiobaccillia = set(ncbi.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 1807140).values())
        joined = ecoli.union(homo).union(acidithiobaccillia)
        expanded = ncbi.expand_sample_taxonomy(sample_set)
        self.assertSetEqual(joined, expanded)

        # make sure we don't get human species
        self.assertNotIn(member=9606, container=expanded)

    def testGetChildren(self):
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        # test case is the family hominidae
        id = 9604

        # should have 4 genus children: gorilla (9592), homo (9605), pongo (9599) and pan (9596)
        expected = {9592, 9605, 9599, 9596}
        self.assertSetEqual(ncbi.get_children(id), expected)

    # def testNumberOfChildren(self):
    #     ncbi = td.ncbi_database_handler(DATA_DIR)
    #
    #     # expect these taxids and all ancestors
    #     sample_set = {9605, 562, 1807140}
    #     ecoli = set(td.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 562, ncbi).values())
    #     homo = set(td.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 9605, ncbi).values())
    #     acidithiobaccillia = set(td.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 1807140, ncbi).values())
    #     joined = ecoli.union(homo).union(acidithiobaccillia)
    #     expanded = td.expand_sample_taxonomy(sample_set, ncbi)
    #
    #     # expect two children of Proteobacteria
    #     self.assertEqual(td.number_of_children(1224, ncbi, 'ncbi', expanded), 2)
    #
    #     # expect zero children of homo
    #     self.assertEqual(td.number_of_children(9605, ncbi, 'ncbi', expanded), 0)
    #
    #     # expect that asking for children of root returns np.inf
    #     # this is to save time, because we'll never filter out root, bacteria, euk, or arch
    #     self.assertEqual(td.number_of_children(1, ncbi, 'ncbi', expanded), np.inf)

    def testGetDescendants(self):
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        id = 9604  # great apes
        # all descendants include 4 genuses/geni (9592, 9605, 9599, 9596)
        # and species:
        # gorilla species (499232, 9593)
        # homo species (1425170, 9606)
        # pan species (9597, 9598)
        # and pongo species (9601, 502961, 9600, 2051901, 9603)
        descendants_exp = {9592, 9605, 9599, 9596, 499232, 9593,
                           1425170, 9606, 9597, 9598,
                           9601, 502961, 9600, 2051901, 9603}
        self.assertSetEqual(ncbi.get_descendants(id), descendants_exp)

    def testGetParents(self):
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        testid = 562 # ecoli
        # parent should be 561, the Escherichia genus
        parent_exp = {561}
        self.assertSetEqual(ncbi.get_parents(testid), parent_exp)

    def testGetAncestors(self):
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        testid = 1692040 # acidiferrobacterales
        # ancestors should be 1224, the proteobacteria phylum, and
        # 1236, the gammaproteobacteria class
        parent_exp = {1224, 1236}
        self.assertSetEqual(ncbi.get_ancestors(testid), parent_exp)

    def testConvertTaxidToName(self):
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        # has to be a list or series
        ecoli = [562]
        # returns list
        name = ncbi.convert_taxid_to_name(ecoli)
        self.assertEqual(name, ['Escherichia coli'])

    def testConvertNameToTaxid(self):
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        name= ['Brassicaceae']
        # returns list
        id = ncbi.convert_name_to_taxid(name)
        self.assertEqual(id, [3700])

    def testConvertNameToTaxid_UnknownNames(self):
        # assigns the ncbi taxid for unidentified (32644) if a name is unknown
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        name = ['Random nonsense']
        id = ncbi.convert_name_to_taxid(name)
        self.assertEqual(id, [32644])

    def testConvertNameToTaxid_Root(self):
        ncbi = td.NCBITaxonomyDb(DATA_DIR)
        name = ['root']
        id = ncbi.convert_name_to_taxid(name)
        self.assertEqual(id, [1])

    def testConvertNameToTaxid_UnipeptResults(self):
        # test a whole bunch of unipept name results
        ncbi = td.NCBITaxonomyDb(DATA_DIR)

        unipept_thaliana = os.path.join(DATA_DIR, 'test', 'unipept_a_thaliana_result.csv')
        with open(unipept_thaliana, mode = 'r') as f:
            f.readline() # ditch column header
            names = [elem.strip('\n') for elem in f.readlines()]
            ids = ncbi.convert_name_to_taxid(names)

        # assure no unknown names (i.e., unipept names are reasonably well supported)
        self.assertEqual(sum([id == 32644 for id in ids]), 0)


if __name__=='__main__':
    unittest.main()
