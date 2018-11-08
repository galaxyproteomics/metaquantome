import subprocess
from tests.testutils import testfile, TTEST_SINFO
import unittest
import pandas as pd
import numpy as np
import os


class TestCLI(unittest.TestCase):
    os.chdir('..')
    def testSingleInt(self):
        out = testfile('cli_out.tab')
        command = '''python3 cli.py expand -m tax --pep_colname peptide --outfile ''' + out
        command += ''' -i metaquantome/data/test/simple_int.tab --tax_file metaquantome/data/test/simple_tax.tab '''
        command += '''--tax_colname "lca" --samps '{"A": ["int"]}' '''
        status = subprocess.call(command, shell=True)
        self.assertEqual(status, 0)

        tax_df = pd.read_csv(out, sep='\t')
        self.assertAlmostEqual(tax_df.query("taxon_name == 'Helicobacter pylori'")['int'].values[0],
                               np.log2(100), places=5)

    def testMultipleInt(self):
        exp_out = testfile('cli_mult_out.tab')
        exp_command = '''python3 cli.py expand -m fn --pep_colname peptide --outfile ''' + exp_out
        exp_command += ''' -i metaquantome/data/test/int_ttest.tab --func_file metaquantome/data/test/multiple_func.tab '''
        exp_command += ''' --func_colname cog --ontology cog ''' + " --samps '" + TTEST_SINFO + "'"
        exp_status = subprocess.call(exp_command, shell=True)
        self.assertEqual(exp_status, 0)

        test_out = testfile('cli_mult_test_out.tab')
        test_command = "python3 cli.py stat -m fn --outfile " + test_out + ' --file ' + exp_out
        test_command += ''' --ontology cog ''' + " --samps '" + TTEST_SINFO + "'"
        test_status = subprocess.call(test_command, shell=True)
        self.assertEqual(test_status, 0)

        test_df = pd.read_csv(test_out, sep="\t", index_col='id')
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(test_df['corrected_p']['C'] > 0.05)
        self.assertTrue(test_df['corrected_p'][['N','D']].le(0.05).all())


if __name__ == '__main__':
    unittest.main()
