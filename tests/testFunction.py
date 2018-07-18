import unittest
from metaquant.runner import metaquant
from metaquant import go
from tests.testutils import testfile
import numpy as np
from metaquant.definitions import DATA_DIR
import os


class TestFunction(unittest.TestCase):
    def testSingleInt(self):
        func=testfile('simple_func.tab')
        int=testfile('simple_int.tab')

        go_df = metaquant('fn', sample_names={'s1': ['int']}, int_file=int, pep_colname='peptide', func_colname='go',
                          func_file=func, ontology='go', test=False)
        self.assertEqual(go_df.loc["GO:0022610"]['int'], np.log2(200))
        self.assertEqual(go_df.loc["GO:0008152"]['int'], np.log2(100))

    def testMultipleInt(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')

        go_df = metaquant('fn', sample_names={'s1': ['int1', 'int2', 'int3']}, func_colname='go',
                          int_file=int, func_file=func, ontology='go')
        self.assertEqual(go_df.loc['GO:0008152']['int1'], np.log2(10))
        self.assertEqual(go_df.loc['GO:0022610']['int2'], np.log2(30))
        # missing values (zeros, nans, NA's, etc) are turned into NaN's
        self.assertTrue(np.isnan(go_df.loc['GO:0000003']['int3']))

    def testDA(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')

        go_df = metaquant('fn', sample_names={'s1': ['int1', 'int2', 'int3'],
                                              's2': ['int4', 'int5', 'int6']}, int_file=int, func_file=func,
                          func_colname='go', ontology='go', test=True)

        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(go_df['corrected_p']['GO:0008152'] > 0.05)
        self.assertTrue(go_df['corrected_p'][['GO:0022610','GO:0000003','GO:0032505']].le(0.05).all())

    def testSlimDown(self):
        func=testfile('func_eggnog.tab')
        int=testfile('int_eggnog.tab')
        go_df = metaquant('fn', sample_names={'NS': ['int737NS', 'int852NS', 'int867NS'],
                                              'WS': ['int737WS', 'int852WS', 'int867WS']}, int_file=int, func_file=func,
                          func_colname='go', ontology='go', slim_down=True, test=True, paired=True)

        # test that all go terms are in slim
        # load slim

        go_dag, go_dag_slim = go.go_database_handler(data_dir=os.path.join(DATA_DIR, 'go'), slim_down=True, overwrite=False)

        returned_gos = set(go_df['id'])

        self.assertTrue(returned_gos.issubset(go_dag_slim.keys()))

    def testCog(self):
        func=testfile('multiple_func.tab')
        int=testfile('multiple_int.tab')

        cog_df = metaquant('fn', sample_names={'s1': ['int1', 'int2', 'int3']}, int_file=int, func_file=func,
                           func_colname='cog', ontology='cog')
        self.assertEqual(cog_df.loc["C"]['s1_mean'], np.log2((10+20+70)/3))
        self.assertEqual(cog_df.loc["N"]['int2'], np.log2(30))

    def testCogTTest(self):
        func=testfile('multiple_func.tab')
        int=testfile('int_ttest.tab')

        cog_df = metaquant('fn', sample_names={'s1': ['int1', 'int2', 'int3'],
                                               's2': ['int4', 'int5', 'int6']}, int_file=int, func_file=func,
                           func_colname='cog', ontology='cog', test=True)
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(cog_df['corrected_p']['C'] > 0.05)
        self.assertTrue(cog_df['corrected_p'][['N','D']].le(0.05).all())


if __name__ == '__main__':
    unittest.main()