import unittest
import os
import pandas as pd
import numpy as np

import metaquantome.databases.NCBITaxonomyDb as td
from metaquantome.util.utils import DATA_DIR, TEST_DIR
from metaquantome.util.testutils import testfile


class TestTaxonomyDatabase(unittest.TestCase):
    ncbi = td.NCBITaxonomyDb(TEST_DIR)

    def testMapIdToDesiredRanks(self):
        # human genus
        # we don't want to get species back, because that's lower in the hierarchy
        human_ancestors = self.ncbi.map_id_to_desired_ranks({'order', 'species', 'genus'}, 9605)
        self.assertDictEqual({'order': 9443, 'genus': 9605}, human_ancestors)

    def testMapIdToDesiredRanks_Unipept(self):
        # take all taxids from unipept sample 7
        # just make sure that there aren't any ranks returned by unipept that are
        # not contained in FULL_TAXONOMY_TREE
        tax = testfile('unipept_sample7_taxonomy.tab')
        df = pd.read_csv(tax, sep='\t')
        num_ids = self.ncbi.convert_name_to_taxid(df['lca'])
        for i in num_ids:
            self.ncbi.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, i)

    def testExpandSampleTaxonomy(self):
        # test that expanded sample set contains all desired ancestors
        # expect these taxids and all ancestors within BASIC_TAXONOMY_TREE
        sample_set = {9605, 562, 1807140}
        ecoli = set(self.ncbi.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 562).values())
        homo = set(self.ncbi.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 9605).values())
        acidithiobaccillia = set(self.ncbi.map_id_to_desired_ranks(td.BASIC_TAXONOMY_TREE, 1807140).values())
        joined = ecoli.union(homo).union(acidithiobaccillia)
        expanded = self.ncbi.expand_sample_taxonomy(sample_set)
        self.assertSetEqual(joined, expanded)
        # make sure we don't get human species
        self.assertNotIn(member=9606, container=expanded)

    def testGetRank(self):
        testid = 2
        rank = self.ncbi.get_rank(testid)
        self.assertEqual(rank, 'superkingdom')

    def testGetChildren(self):
        # test case is the family hominidae
        id = 9604
        # should have 4 genus children: gorilla (9592), homo (9605), pongo (9599) and pan (9596)
        expected = {9592, 9605, 9599, 9596}
        self.assertSetEqual(self.ncbi.get_children(id), expected)

    def testGetDescendants(self):
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
        self.assertSetEqual(self.ncbi.get_descendants(id), descendants_exp)

    def testGetParents(self):
        testid = 562 # ecoli
        # parent should be 561, the Escherichia genus
        parent_exp = {561}
        self.assertSetEqual(self.ncbi.get_parents(testid), parent_exp)

    def testGetAncestors(self):
        testid = 1692040 # acidiferrobacterales
        # ancestors should be 1224, the proteobacteria phylum, and
        # 1236, the gammaproteobacteria class
        parent_exp = {1224, 1236}
        self.assertSetEqual(self.ncbi.get_ancestors(testid), parent_exp)

    def testConvertTaxidToName(self):
        # has to be a list or series
        ecoli = [562]
        # returns list
        name = self.ncbi.convert_taxid_to_name(ecoli)
        self.assertEqual(name, ['Escherichia coli'])

    def testConvertNameToTaxid(self):
        name= ['Brassicaceae']
        # returns list
        id = self.ncbi.convert_name_to_taxid(name)
        self.assertEqual(id, [3700])

    def testConvertNameToTaxid_UnknownNames(self):
        # assigns the ncbi taxid for unidentified (32644) if a name is unknown
        name = ['Random nonsense']
        id = self.ncbi.convert_name_to_taxid(name)
        self.assertEqual(id, [np.nan])

    def testConvertNameToTaxid_Root(self):
        name = ['root']
        id = self.ncbi.convert_name_to_taxid(name)
        self.assertEqual(id, [1])

    def testConvertNameToTaxid_UnipeptResults(self):
        # test a whole bunch of unipept name results
        unipept_thaliana = os.path.join(DATA_DIR, 'test', 'unipept_a_thaliana_result.csv')
        with open(unipept_thaliana, mode='r') as f:
            f.readline() # ditch column header
            names = [elem.strip('\n') for elem in f.readlines()]
            ids = self.ncbi.convert_name_to_taxid(names)
        # assure no unknown names (i.e., unipept names are reasonably well supported)
        self.assertEqual(sum([id == np.nan for id in ids]), 0)


if __name__=='__main__':
    unittest.main()
