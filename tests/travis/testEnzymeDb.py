import unittest

from metaquantome.databases.EnzymeDb import EnzymeDb
from metaquantome.util.utils import TEST_DIR


class TestEC(unittest.TestCase):
    ec = EnzymeDb(TEST_DIR)

    def testAnnotateEcDatabase(self):
        testdb = {'1.-.-.-': "a description",
                  '1.2.-.-': "another description",
                  '1.2.3.-': "more specific enzyme",
                  '1.2.3.4': "most specific enzyme"}
        newdb = self.ec._annotate_enzyme_db(testdb)
        expected0 = {'descript': testdb['1.-.-.-'],
                     'depth': 0,
                     'levels': ['1', '-', '-', '-'],
                     'id': '1.-.-.-'}
        expected1 = {'descript': testdb['1.2.-.-'],
                     'depth': 1,
                     'levels': ['1', '2', '-', '-'],
                     'id': '1.2.-.-'}
        expected2 = {'descript': testdb['1.2.3.-'],
                     'depth': 2,
                     'levels': ['1', '2', '3', '-'],
                     'id': '1.2.3.-'}
        expected3 = {'descript': testdb['1.2.3.4'],
                     'depth': 3,
                     'levels': ['1', '2', '3', '4'],
                     'id': '1.2.3.4'}
        all_expected = {'1.-.-.-': expected0,
                        '1.2.-.-': expected1,
                        '1.2.3.-': expected2,
                        '1.2.3.4': expected3}
        self.assertDictEqual(newdb, all_expected)

    def testAssignDepth(self):
        test_ecid0 = '1.-.-.-'
        test_ecid2 = '4.5.2.-'
        self.assertEqual(self.ec._assign_depth(test_ecid0), 0)
        self.assertEqual(self.ec._assign_depth(test_ecid2), 2)

    def testAssignLevels(self):
        test_ecid0 = '1.-.-.-'
        test_ecid2 = '4.5.2.-'
        self.assertEqual(self.ec._split_ec(test_ecid0), ['1', '-', '-', '-'])
        self.assertEqual(self.ec._split_ec(test_ecid2), ['4', '5', '2', '-'])

    def testGetChildren(self):
        # Oxidoreductases, Acting on the CH-OH group of donors, With a disulfide as acceptor
        test_id = '1.1.4.-'
        # expect two children, 1.1.4.1 and 1.1.4.2
        exp_children = {'1.1.4.1', '1.1.4.2'}
        self.assertSetEqual(self.ec.get_children(test_id), exp_children)

    def testGetDescendants(self):
        test_id = '6.5.-.-'  # ligases, forming phosporic ester bonds
        # expect 6.5.1.[1 through 8]
        exp_children = {'6.5.1.' + str(x) for x in range(1, 9)}
        exp_children.update({'6.5.1.-'})  # don't forget about intermediary
        self.assertSetEqual(self.ec.get_descendants(test_id), exp_children)

    def testGetParents(self):
        test_id = '6.5.-.-'  # ligases, forming phosporic ester bonds
        # expect 6.-.-.-
        exp_parent = {'6.-.-.-'}
        self.assertSetEqual(self.ec.get_parents(test_id), exp_parent)

        # if no parents, returns empty set
        test_id2 = '6.-.-.-'
        exp_parent2 = set()
        self.assertSetEqual(self.ec.get_parents(test_id2), exp_parent2)

    def testGetAncestors(self):
        test_id1 = '6.5.1.-'  # ligases, forming phosporic ester bonds
        # expect 6.-.-.-
        exp_anc1 = {'6.5.-.-', '6.-.-.-'}
        self.assertSetEqual(self.ec.get_ancestors(test_id1), exp_anc1)

        # if no ancestors, returns empty set
        test_id2 = '6.-.-.-'
        exp_parent2 = set()
        self.assertSetEqual(self.ec.get_ancestors(test_id2), exp_parent2)

        # includes root
        test_id3 = '1.1.4.1'
        exp_anc3 = {'1.-.-.-', '1.1.-.-', '1.1.4.-'}
        self.assertSetEqual(self.ec.get_ancestors(test_id3), exp_anc3)


if __name__ == '__main__':
    unittest.main()
