import unittest

from metaquantome.databases import GeneOntologyDb as godb
from metaquantome.util.utils import TEST_DIR


class TestGeneOntologyDb(unittest.TestCase):
    db = godb.GeneOntologyDb(TEST_DIR, slim_down=True)

    def testMapIdToSlim(self):
        # the case where the test term has a parent in the slim set
        testid = "GO:0009343" # biotin carboxylase complex
        # closest ancestor in slim is GO:1902494, catalytic complex
        exp_closest = "GO:1902494"
        obs_closest = self.db.map_id_to_slim(testid)
        self.assertEqual(exp_closest, obs_closest)

        # case where test term has a grandparent in the slim set
        test_grand_id = "GO:0032284"  # GO:0032284 plastid biotin carboxylase complex
        # closest ancestor in slim is again GO:1902494, catalytic complex
        exp_grand_closest = "GO:1902494"
        obs_grand_closest = self.db.map_id_to_slim(test_grand_id)
        self.assertEqual(exp_grand_closest, obs_grand_closest)

        # case where go term is not in full set
        test_gibberish = "gibberish"
        exp_result = 'unknown'
        obs_result = self.db.map_id_to_slim(test_gibberish)
        self.assertEqual(exp_result, obs_result)

    def testMapSetToSlim(self):
        # use the same terms as in testMapIdToSlim
        test_set = {"GO:0009343", "GO:0032284", "gibberish"}
        exp_result = {"GO:0009343": "GO:1902494",
                      "GO:0032284": "GO:1902494",
                      "gibberish": 'unknown'}
        obs_result = self.db.map_set_to_slim(test_set)
        self.assertDictEqual(exp_result, obs_result)

    def testGetChildren(self):
        testid = "GO:0098632"  # cell-cell adhesion mediator activity
        # expected children
        exp_children = {"GO:0086080",
                        "GO:0098641"}
        self.assertSetEqual(exp_children, self.db.get_children(testid))

    def testGetChildren_Unknown(self):
        testid = "notagoterm"
        self.assertSetEqual(self.db.get_children(testid), set())

    def testGetDescendants(self):
        testid = "GO:0007044"  # cell-substrate junction assembly
        # two children, three descendants:
        # GO:0007045
        #   - GO:0048041
        # GO:0031581
        exp_descendants = {"GO:0007045",
                           "GO:0048041",
                           "GO:0031581"}
        act_descendants = self.db.get_descendants(testid)
        self.assertSetEqual(exp_descendants, act_descendants)

    def testGetDescendants_Unknown(self):
        testid = "notagoterm"
        self.assertSetEqual(self.db.get_descendants(testid), set())

    def testGetParents(self):
        testid = "GO:0044406"  # cell-cell adhesion mediator activity
        # expected parents
        exp_parents = {"GO:0022610",
                       "GO:0051704"}
        self.assertSetEqual(exp_parents, self.db.get_parents(testid))

    def testGetParents_Unknown(self):
        testid = "notagoterm"
        self.assertSetEqual(self.db.get_parents(testid), set())

    def testGetAncestors(self):
        testid = "GO:0016043"  # cellular component organization
        # two parents and one grandparent
        exp_ancestors = {"GO:0009987",
                         "GO:0071840",
                         "GO:0008150"}
        act_ancestors = self.db.get_ancestors(testid)
        self.assertSetEqual(exp_ancestors, act_ancestors)

    def testGetAncestors_Unknown(self):
        testid = "notagoterm"
        self.assertSetEqual(self.db.get_ancestors(testid), set())


if __name__=="__main__":
    unittest.main()