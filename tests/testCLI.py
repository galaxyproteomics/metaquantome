import subprocess
from tests.testutils import testfile
import unittest
import os
import pandas as pd
import numpy as np
from metaquant.__main__ import main


class TestCLI(unittest.TestCase):
    def testSingleInt(self):
        out=testfile('cli_out.tab')

        os.chdir("..")

        command = '''python3 metaquant/__main__.py -m tax --pep_colname peptide --outfile ''' + out
        command += ''' -i metaquant/data/test/simple_int.tab --tax_file metaquant/data/test/simple_tax.tab '''
        command += '''--tax_colname "lca" --samps '{"A": ["int"]}' '''
        subprocess.call(command, shell=True)

        tax_df = pd.read_csv(out, sep='\t')

        self.assertAlmostEqual(tax_df.query("taxon_name == 'Helicobacter pylori'")['int'].values[0], np.log2(100), places=5)
