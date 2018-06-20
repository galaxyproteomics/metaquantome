import unittest
from src import io
import os
from definitions import DATA_DIR


class TestIO(unittest.TestCase):
    def testTaxonomyIn(self):
        taxin = os.path.join(DATA_DIR, 'test', 'unipept_a_thaliana_result_w_pep.tab')
        df = io.read_taxonomy_table(taxin, 'peptide', 'lca')
        # the first peptide is assigned to 'root', which is taxid 1
        # is translated to taxids
        self.assertEqual(df['lca'][0], 1)

    def testFunctionIn(self):
        funcin = os.path.join(DATA_DIR, 'test', 'multiple_func.tab')
        df = io.read_function_table(funcin, 'peptide', 'go')
        self.assertEqual(df['go']['A'], 'GO:0008152')

    def testMerge(self):
        taxin = os.path.join(DATA_DIR, 'test', 'tax_join.tab')
        funcin = os.path.join(DATA_DIR, 'test', 'func_join.tab')
        int_in = os.path.join(DATA_DIR, 'test', 'int_join.tab')

        samp_grps = io.SampleGroups({'1':['int1', 'int2']})

        dfs_joined = io.read_and_join_files('taxfn', pep_colname='peptide',
                                            int_file=int_in, samp_groups=samp_grps,
                                            tax_file=taxin, func_file=funcin, tax_colname='lca', func_colname='go')

        self.assertSetEqual(set(dfs_joined), {'int1', 'int2', 'lca', 'go'})

    def testSampInfo(self):
        # test that info file is correctly parsed
        sinfo = os.path.join(DATA_DIR, 'test', 'samp_info.tab')
        sample_groups = io.read_samp_info(sinfo)
        self.assertDictEqual(sample_groups, {'A': ['A1', 'A2'], 'B': ['B1', 'B2']})

        # test that we can create SampleGroups() with this
        samp_grps = io.SampleGroups(sample_groups)
        self.assertEqual(samp_grps.grp_names, ['A', 'B'])

        # note that the order of all_intcols is random (so we compare sets)
        self.assertEqual(set(samp_grps.all_intcols), {'A1', 'A2', 'B1', 'B2'})

