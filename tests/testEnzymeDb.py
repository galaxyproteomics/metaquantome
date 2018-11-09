import unittest
from metaquantome.databases.EnzymeDb import EnzymeDb
from metaquantome.util.utils import DATA_DIR
import os
import shutil


class TestEC(unittest.TestCase):
    TEST_DIR = os.path.join(DATA_DIR, 'test', 'ec_cache')  # downloaded 8/28/18
    ec = EnzymeDb(TEST_DIR)

    def testEnzymeDatabaseHandler(self):
        # make tmp dir for testing download
        tmp_dir = os.path.join(DATA_DIR, 'tmp_test_data_dwnld')
        os.mkdir(tmp_dir)
        try:
            enzyme_db = EnzymeDb(tmp_dir, False)
            expected_contents = [os.path.join(tmp_dir, file)
                                 for file in ['enzclass.txt', 'enzyme.dat', 'ec_id.json', 'enzclass.json']]
            for content in expected_contents:
                self.assertTrue(os.path.exists(content))
            # make sure parsed correctly
            # this is from enzyme.dat
            self.assertEqual(enzyme_db.ecdb['1.2.3.4']['descript'], 'Oxalate oxidase.')

            # from enzclass.txt
            self.assertEqual(enzyme_db.ecdb['6.1.-.-']['descript'], 'Forming carbon-oxygen bonds.')
        finally:
            shutil.rmtree(tmp_dir)

    def testAnnotateEcDatabase(self):
        testdb = {'1.-.-.-': "a description",
                  '1.2.-.-': "another description",
                  '1.2.3.-': "more specific enzyme",
                  '1.2.3.4': "most specific enzyme"}
        newdb = self.ec._annotate_ec_database(testdb)
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
        self.assertEqual(self.ec._assign_levels(test_ecid0), ['1', '-', '-', '-'])
        self.assertEqual(self.ec._assign_levels(test_ecid2), ['4', '5', '2', '-'])

    # def testExpandEC(self):
    #     test_ec = '1.2.-.-'
    #     df = ec.expand_ec(test_ec)
    #     # should have two levels (the first two)
    #     self.assertEqual(df.loc[ec.LEVEL_NAMES[0]], '1.-.-.-')
    #     self.assertEqual(df.loc[ec.LEVEL_NAMES[1]], '1.2.-.-')
    #     self.assertTrue(ec.LEVEL_NAMES[2] not in df.index)

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
        test_id = '6.5.1.-'  # ligases, forming phosporic ester bonds
        # expect 6.-.-.-
        exp_parent = {'6.5.-.-', '6.-.-.-'}
        self.assertSetEqual(self.ec.get_ancestors(test_id), exp_parent)

        # if no ancestors, returns empty set
        test_id2 = '6.-.-.-'
        exp_parent2 = set()
        self.assertSetEqual(self.ec.get_ancestors(test_id2), exp_parent2)
