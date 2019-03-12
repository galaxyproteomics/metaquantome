import unittest
import pandas as pd

from metaquantome.classes.AnnotationNode import AnnotationNode
from metaquantome.classes.AnnotationHierarchy import AnnotationHierarchy
from metaquantome.databases.NCBITaxonomyDb import NCBITaxonomyDb
from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
from metaquantome.databases.EnzymeDb import EnzymeDb
from metaquantome.util.utils import TEST_DIR


class TestAnnotationHierarchyNcbi(unittest.TestCase):

    def _create_sapiens_db(self):
        db = NCBITaxonomyDb(TEST_DIR)
        sample_set = {9604, 9605, 9606}  # hominidae (family), homo (genus), homo sapiens (species)
        ah = AnnotationHierarchy(db, sample_set, 'samp1')
        return ah, sample_set

    def testInit(self):
        ah, sample_set = self._create_sapiens_db()
        self.assertIsInstance(ah.db, NCBITaxonomyDb)
        self.assertSetEqual(ah.sample_set, sample_set)
        self.assertDictEqual(ah.nodes, dict())

    def testUpdateNode(self):
        ah, sample_set = self._create_sapiens_db()
        # one sample child
        testid = 9605
        intensity = 200
        ah._add_node_with_ancestors(testid, intensity)

        testid2 = 9606
        ah._add_node_with_ancestors(testid2, intensity)

        ah._define_sample_children()
        updated_node = ah.nodes[testid]
        self.assertIsInstance(updated_node, AnnotationNode)
        self.assertEqual(updated_node.intensity, intensity*2)
        self.assertEqual(updated_node.n_sample_children, 1)

    def testAggregateNodes(self):
        ah, sample_set = self._create_sapiens_db()
        testids = [9604, 9605, 9606]
        test_intensities = [500, 200, 300]
        for i in range(0, 3):
            ah._add_node_with_ancestors(testids[i], test_intensities[i])
        self.assertEqual(ah.nodes[9604].intensity, 1000)


class TestAnnotationHierarchyGO(unittest.TestCase):
    def _create_go_db(self):
        db = GeneOntologyDb(TEST_DIR)
        sample_set = {'GO:0008150',  # biological process
                      'GO:0008283',  # cell proliferation (child of BP)
                      'GO:0033687',  # osteoblast proliferation (child of cell pro)
                      'GO:0036093',  # germ cell proliferation (child of cell pro)
                      'GO:0022414',  # reproductive process (child of BP)
                      'GO:1903046',  # meiotic cell cycle process (child of rep pro)
                      'GO:0051026'}  # chiasma assembly, child of meiotic
        ah = AnnotationHierarchy(db, sample_set, 'samp1')
        return ah, sample_set

    def testInit(self):
        ah, sample_set = self._create_go_db()
        self.assertIsInstance(ah.db, GeneOntologyDb)
        self.assertSetEqual(ah.sample_set, sample_set)
        self.assertDictEqual(ah.nodes, dict())

    def testUpdateNode(self):
        ah, sample_set = self._create_go_db()
        # one sample child
        testid = 'GO:0051026'
        intensity = 100
        ah._add_node_with_ancestors(testid, intensity)
        updated_node = ah.nodes[testid]
        self.assertIsInstance(updated_node, AnnotationNode)
        self.assertEqual(updated_node.intensity, intensity)
        ah._define_sample_children()
        self.assertEqual(updated_node.n_sample_children, 0)

    def testAggregateNodes(self):
        ah, sample_set = self._create_go_db()
        testids = ['GO:0008150',  # biological process
                   'GO:0008283',  # cell proliferation (child of BP)
                   'GO:0033687',  # osteoblast proliferation (child of cell pro)
                   'GO:0036093',  # germ cell proliferation (child of cell pro and rep pro)
                   'GO:0022414',  # reproductive process (child of BP)
                   'GO:1903046',  # meiotic cell cycle process (child of rep pro)
                   'GO:0051026']  # chiasma assembly, child of meiotic
        test_intensities = [0, 0, 0, 100, 50, 200, 300]
        for i in range(0, len(test_intensities)):
            ah._add_node_with_ancestors(testids[i], test_intensities[i])
        self.assertEqual(ah.nodes['GO:0022414'].intensity, 650)


class TestAnnotationHierarchyEc(unittest.TestCase):
    def _create_ec_db(self):
        db = EnzymeDb(TEST_DIR)
        sample_set = {'1.1.4.-',
                      '1.1.4.1',
                      '1.1.4.2',
                      '6.5.-.-',
                      '6.-.-.-'}
        ah = AnnotationHierarchy(db, sample_set, 'samp1')
        return ah, sample_set

    def testInit(self):
        ah, sample_set = self._create_ec_db()
        self.assertIsInstance(ah.db, EnzymeDb)
        self.assertSetEqual(ah.sample_set, sample_set)
        self.assertDictEqual(ah.nodes, dict())

    def testUpdateNode(self):
        ah, sample_set = self._create_ec_db()
        # one sample child
        testids = ['1.1.4.-', '1.1.4.1', '1.1.4.2']
        intensity = 100
        for i in testids:
            ah._add_node_with_ancestors(i, intensity)
        updated_node = ah.nodes[testids[0]]
        print(ah.nodes['1.1.-.-'].intensity)
        ah._define_sample_children()
        self.assertIsInstance(updated_node, AnnotationNode)
        self.assertEqual(updated_node.intensity, intensity*3)
        self.assertEqual(updated_node.n_sample_children, 2)

    def testAggregateNodes(self):
        ah, sample_set = self._create_ec_db()
        testids = ['1.1.4.-',
                   '1.1.4.1',
                   '1.1.4.2',
                   '6.5.-.-']
        test_intensities = [500, 200, 300, 0]
        for i in range(0, 3):
            ah._add_node_with_ancestors(testids[i], test_intensities[i])
        self.assertEqual(ah.nodes['1.1.4.-'].intensity, 1000)

    def testToDataframe(self):
        ah, sample_set = self._create_ec_db()
        test_set = ['1.1.4.-',
                    '1.1.4.1',
                    '1.1.4.2',
                    '1.1.4.-',
                    '1.1.4.1',
                    '1.1.4.2',
                    '6.5.-.-',
                    '6.5.-.-',
                    '6.-.-.-',
                    '6.-.-.-']
        # set to one, so it's equal to number of peptides
        test_intensity = 1
        for id in test_set:
            ah._add_node_with_ancestors(id, test_intensity)

        # expanded sample set is all nodes
        ah._define_sample_children()

        # the sample set is as below:
        #         sample_set = {'1.1.4.-',
        #                       '1.1.4.1',
        #                       '1.1.4.2',
        #                       '6.5.-.-',
        #                       '6.-.-.-'}
        # ah.get_informative_nodes(0, 0)
        # expected
        exp_df = pd.DataFrame({'samp1': [6, 6, 6, 2, 2, 4, 2],
                               'samp1_n_peptide': [6, 6, 6, 2, 2, 4, 2],
                               'samp1_n_samp_children': [1, 1, 2, 0, 0, 1, 0]},
                               index= ['1.-.-.-',
                                       '1.1.-.-',
                                       '1.1.4.-',
                                       '1.1.4.1',
                                       '1.1.4.2',
                                       '6.-.-.-',
                                       '6.5.-.-']).sort_index()
        df = ah.to_dataframe().sort_index()
        self.assertTrue(df.equals(exp_df))


if __name__ == '__main__':
    unittest.main()
