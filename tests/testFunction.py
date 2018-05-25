import unittest
import metaquant
import os
from definitions import DATA_DIR
from src import go


class TestFunction(unittest.TestCase):
    def testSingleInt(self):
        datfile=os.path.join(DATA_DIR, 'test', 'function_single_int.tab')
        go_df = metaquant.metaquant('fn', datfile,
                                    func_colname='go', sample_names={'s1': ['int']}, test=False,
                                    ontology="GO")
        self.assertEqual(go_df.loc["GO:0022610"]['int'], 2/3)
        self.assertEqual(go_df.loc["GO:0008152"]['int'], 1 / 3)

    def testMultipleInt(self):
        datfile=os.path.join(DATA_DIR, 'test', 'function_multiple_int.tab')
        go_df = metaquant.metaquant('fn', datfile, func_colname='go',
                                    sample_names={'s1': ['int1', 'int2', 'int3']}, test=False,
                                    ontology="GO")
        self.assertEqual(go_df.loc['GO:0008152']['int1'], 1/3)
        self.assertEqual(go_df.loc['GO:0022610']['int2'], 3/5)
        self.assertEqual(go_df.loc['GO:0008152']['int3'], 0.7)

    def testRedundant(self):
        """
        test that if a parent and child term are both present, the parent doesn't get double
        """
        datfile=os.path.join(DATA_DIR, 'test', 'function_multiple_go.tab')
        godf = metaquant.metaquant('fn', datfile,
                                   func_colname='go',
                                   sample_names={'s1': ['int1', 'int2', 'int3']},
                                   test=False,
                                   ontology="GO")
        self.assertEqual(godf.loc['GO:0008152']['int1'], 1/3)
        self.assertEqual(godf.loc['GO:0022610']['int2'], 3/5)
        self.assertEqual(godf.loc['GO:0008152']['int3'], 0.7)

    def testEggnogOutput(self):
        datfile=os.path.join(DATA_DIR, 'test', 'function_eggnog_gos.tabular')
        go_df = metaquant.metaquant('fn', datfile, func_colname='go',
                                    ontology="GO",
                                    sample_names={'s1': ['int737WS', 'int737NS', 'int852WS',
                                                         'int852NS', 'int867WS', 'int867NS']},
                                    test=False)
        self.assertTrue(go_df['int852WS'].le(1).all())

    def testDA(self):
        datfile = os.path.join(DATA_DIR, 'test', 'function_multiple_int_ttests.tab')
        go_df = metaquant.metaquant('fn', datfile,
                                    func_colname='go',
                                    ontology="GO",
                                    sample_names={'s1': ['int1', 'int2', 'int3'],
                                                  's2': ['int4', 'int5', 'int6']},
                                    test=True)
        # make sure all are less than 0.05
        self.assertTrue(go_df['corrected_p'].le(0.05).all())

    def testSlimDown(self):
        datfile=os.path.join(DATA_DIR, 'test', 'function_eggnog_gos.tabular')
        go_df = metaquant.metaquant('fn', datfile, func_colname='go',
                                    ontology="GO",
                                    sample_names={'NS': ['int737NS', 'int852NS', 'int867NS'],
                                                  'WS': ['int737WS', 'int852WS', 'int867WS']},
                                    test=True, slim_down=True,
                                    paired=True)

        # test that all go terms are in slim
        # load slim
        go_dag, go_dag_slim = go.load_obos(slim_down=True)

        returned_gos = set(go_df['id'])

        self.assertTrue(returned_gos.issubset(go_dag_slim.keys()))

    def testCog(self):
        datfile=os.path.join(DATA_DIR, 'test', 'function_cog_single_int.tab')
        cog_df = metaquant.metaquant('fn', datfile,
                                    func_colname="cog",
                                    ontology="cog",
                                    sample_names={'s1' : ['int']},
                                    test=False)
        self.assertEqual(cog_df.loc["O"]['int'], 5/7)
        self.assertEqual(cog_df.loc["C"]['int'], 2/7)

    def testCogTTest(self):
        datfile=os.path.join(DATA_DIR, 'test', 'function_cog_multiple_int_ttest.tab')
        cog_df = metaquant.metaquant('fn', datfile,
                                    func_colname="cog",
                                    ontology="cog",
                                    sample_names={'samp1': ['int1', 'int2', 'int3'],
                                                  'samp2': ['int4', 'int5', 'int6']},
                                    test=True)
        # make sure all are less than 0.05
        self.assertTrue(cog_df['corrected_p'].le(0.05).all())


    def testCogFiltering(self):
        # the cog category O has only 1 observation in sample1, so it should be filtered out at threshold 2
        datfile=os.path.join(DATA_DIR, 'test', 'function_cog_multiple_int_ttest_filtering.tab')
        cog_df = metaquant.metaquant('fn', datfile,
                                    func_colname="cog",
                                    ontology="cog",
                                    sample_names={'samp1': ['int1', 'int2', 'int3'],
                                                  'samp2': ['int4', 'int5', 'int6']},
                                    test=True,
                                    threshold=2)
        self.assertTrue("O" not in set(cog_df.index) and "C" in set(cog_df.index))


if __name__ == '__main__':
    unittest.main()