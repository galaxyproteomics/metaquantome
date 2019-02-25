import unittest
import os

import metaquantome.classes.SampleGroups
import metaquantome.util.expand_io
from metaquantome.util.utils import DATA_DIR
from metaquantome.util.testutils import testfile


class TestIO(unittest.TestCase):

    def testTaxonomyIn(self):
        taxin = os.path.join(DATA_DIR, 'test', 'unipept_a_thaliana_result_w_pep.tab')
        df = metaquantome.util.expand_io.read_taxonomy_table(taxin, 'peptide', 'lca')
        # the first peptide is assigned to 'root'
        self.assertEqual(df['lca'][0], 'root')

    def testFunctionIn(self):
        funcin = os.path.join(DATA_DIR, 'test', 'multiple_func.tab')
        df = metaquantome.util.expand_io.read_function_table(funcin, 'peptide', func_colname='go')
        self.assertEqual(df['go']['A'], 'GO:0008152')

    def testMerge(self):
        taxin = os.path.join(DATA_DIR, 'test', 'tax_join.tab')
        funcin = os.path.join(DATA_DIR, 'test', 'func_join.tab')
        int_in = os.path.join(DATA_DIR, 'test', 'int_join.tab')
        sinfo = '{"1":["int1", "int2"]}'
        samp_grps = metaquantome.classes.SampleGroups.SampleGroups(sinfo)
        dfs_joined = metaquantome.util.expand_io.read_and_join_files('ft', pep_colname_int='peptide',
                                                                     pep_colname_func=None, pep_colname_tax=None,
                                                                     samp_grps=samp_grps, int_file=int_in,
                                                                     tax_file=taxin, func_file=funcin,
                                                                     func_colname='go', tax_colname='lca')

        self.assertSetEqual(set(dfs_joined), {'int1', 'int2', 'lca', 'go'})

    def testNopepIn(self):
        nopep = testfile('nopep.tab')
        sinfo = '{"s1":["int1", "int2", "int3"]}'
        samp_grps = metaquantome.classes.SampleGroups.SampleGroups(sinfo)

        go_df = metaquantome.util.expand_io.read_nopep_table(mode='f', samp_grps=samp_grps,
                                                             file=nopep, func_colname='cog')
        # test that go is in columns
        self.assertIn('cog', go_df.keys())
        # test that there are only two rows
        # the row with an NA cog was filtered out
        self.assertEqual(go_df.shape[0], 2)

if __name__=='__main__':
    unittest.main()
