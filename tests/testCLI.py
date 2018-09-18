import subprocess
from tests.testutils import testfile
import unittest
import pandas as pd
import numpy as np
import os

class TestCLI(unittest.TestCase):
    def testSingleInt(self):
        out = testfile('cli_out.tab')
        os.chdir('..')
        command = '''python3 metaquant/cli.py -m tax --pep_colname peptide --outfile ''' + out
        command += ''' -i metaquant/data/test/simple_int.tab --tax_file metaquant/data/test/simple_tax.tab '''
        command += '''--tax_colname "lca" --samps '{"A": ["int"]}' --min_peptides 0 --min_children_non_leaf 0'''
        status = subprocess.call(command, shell=True)
        self.assertEqual(status, 0)

        tax_df = pd.read_csv(out, sep='\t')
        self.assertAlmostEqual(tax_df.query("taxon_name == 'Helicobacter pylori'")['int'].values[0],
                               np.log2(100), places=5)


if __name__ == '__main__':
    unittest.main()
