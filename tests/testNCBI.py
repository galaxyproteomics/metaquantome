import unittest
from metaquant import phylo_tree
import os
from metaquant.definitions import DATA_DIR


class TestNCBI(unittest.TestCase):
    def testGoUpTree(self):
        ncbi = phylo_tree.load_ncbi()
        # human
        human_ancestors = phylo_tree.get_desired_ranks_from_lineage({'order', 'species', 'genus'}, 9606, ncbi)
        self.assertDictEqual({'order': 9443, 'species' : 9606, 'genus': 9605}, human_ancestors)

    def testTranslator(self):
        ncbi = phylo_tree.load_ncbi()
        # has to be a list or series
        ecoli = [562]

        # returns list
        name = phylo_tree.convert_taxid_to_name(ecoli, ncbi)
        self.assertEqual(name, ['Escherichia coli'])

    def testTranslateName2Id(self):
        ncbi = phylo_tree.load_ncbi()

        name= ['Brassicaceae']

        # returns list
        id = phylo_tree.convert_name_to_taxid(name, ncbi)
        self.assertEqual(id, [3700])

    def testUnknownNames(self):
        # assigns the ncbi taxid for unassigned (32644) if a name is unknown

        ncbi = phylo_tree.load_ncbi()

        name = ['Random nonsense']

        id = phylo_tree.convert_name_to_taxid(name, ncbi)
        self.assertEqual(id, [12908])


    def testRoot(self):
        ncbi = phylo_tree.load_ncbi()

        name = ['root']

        id = phylo_tree.convert_name_to_taxid(name, ncbi)
        self.assertEqual(id, [1])

    def testUnipeptResults(self):
        # test a whole bunch of unipept name results
        ncbi = phylo_tree.load_ncbi()

        unipept_thaliana = os.path.join(DATA_DIR, 'test', 'unipept_a_thaliana_result.csv')
        with open(unipept_thaliana, mode = 'r') as f:
            f.readline() # ditch column header
            names = [elem.strip('\n') for elem in f.readlines()]
            ids = phylo_tree.convert_name_to_taxid(names, ncbi)

        # assure no unknown names (i.e., unipept names are reasonably well supported)
        self.assertEqual(sum([id == 32644 for id in ids]), 0)



if __name__=='__main__':
    unittest.main()
