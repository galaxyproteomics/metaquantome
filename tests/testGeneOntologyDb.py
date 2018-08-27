import unittest
from metaquant import GeneOntologyDb as godb
from metaquant.definitions import DATA_DIR
import os
import shutil


class TestGeneOntologyDb(unittest.TestCase):
    TEST_DIR = os.path.join(DATA_DIR, 'test', 'go_cache')  # downloaded 8/27/18
    db = godb.GeneOntologyDb(TEST_DIR)

    def testDownloadGO(self):
        tmp_dir = os.path.join(DATA_DIR, 'tmp_go_data_dwnld')
        os.mkdir(tmp_dir)
        try:
            go = godb.GeneOntologyDb(tmp_dir)
            expected_contents = [os.path.join(tmp_dir, file)
                                 for file in ['go-basic.obo', 'goslim_generic.obo']]
            for content in expected_contents:
                self.assertTrue(os.path.exists(content))
            # make sure parsed correctly
            self.assertEqual('biological_process', go.gofull['GO:0008150'].name)
        finally:
            shutil.rmtree(tmp_dir)

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
