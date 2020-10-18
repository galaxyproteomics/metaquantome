import subprocess
import unittest
import pandas as pd
import numpy as np
import os

from metaquantome.util.utils import TEST_DIR
from metaquantome.util.testutils import testfile, TTEST_SINFO


class TestCLI(unittest.TestCase):
    def testSingleInt(self):
        out = testfile('cli_out.tab')
        command = '''python3 metaquantome/cli.py expand -m t --pep_colname_int peptide --pep_colname_tax peptide ''' +\
            '''--outfile ''' + out
        command += ''' -i metaquantome/data/test/simple_int.tab --tax_file metaquantome/data/test/simple_tax.tab '''
        command += '''--tax_colname "lca" --samps '{"A": ["int"]}' '''
        command += '''--data_dir ''' + TEST_DIR
        status = subprocess.call(command, shell=True)
        self.assertEqual(status, 0)

        tax_df = pd.read_csv(out, sep='\t')
        self.assertAlmostEqual(tax_df.query("taxon_name == 'Helicobacter pylori'")['int'].values[0],
                               np.log2(100), places=5)

    def testMultipleInt(self):
        exp_out = testfile('cli_mult_out.tab')
        exp_command = '''python3 metaquantome/cli.py expand -m f --pep_colname_int peptide --pep_colname_func peptide ''' +\
            '''--outfile ''' + exp_out
        exp_command += ''' -i metaquantome/data/test/int_ttest.tab --func_file metaquantome/data/test/multiple_func.tab '''
        exp_command += ''' --func_colname cog --ontology cog ''' + " --samps '" + TTEST_SINFO + "' "
        exp_command += '''--data_dir ''' + TEST_DIR
        exp_status = subprocess.call(exp_command, shell=True)
        self.assertEqual(exp_status, 0)

        test_out = testfile('cli_mult_test_out.tab')
        test_command = "python3 metaquantome/cli.py stat -m f --outfile " + test_out + ' --file ' + exp_out + ' --control_group s2 ' 
        test_command += ''' --ontology cog ''' + " --samps '" + TTEST_SINFO + "'" + ' --parametric True '
        test_status = subprocess.call(test_command, shell=True)
        self.assertEqual(test_status, 0)

        test_df = pd.read_csv(test_out, sep="\t", index_col='id')
        # make sure false is > 0.05 and trues are less than 0.05
        self.assertTrue(test_df['corrected_p_s1_over_s2']['C'] > 0.05)
        self.assertTrue(test_df['corrected_p_s1_over_s2'][['N','D']].le(0.05).all())

    def testViz(self):
        infile = testfile('taxonomy_write_simple.tab')
        imgfile = testfile('cli_bar_viz.png')
        cmd = ' '.join([
            'python3 metaquantome/cli.py viz -m t --plottype bar --infile',
            infile,
            '--img',
            imgfile,
            """--samps '{"samp1": ["int"]}'""",
            '--nterms 2 --meancol samp1_mean --target_rank genus'
        ])
        test_status = subprocess.call(cmd, shell=True)
        self.assertEqual(test_status, 0)

    def testVizTabfile(self):
        infile = testfile('taxonomy_write_simple.tab')
        imgfile = testfile('cli_bar_viz2.png')
        tabfile = testfile("tmp")
        cmd = ' '.join([
            'python3 metaquantome/cli.py viz -m t --plottype bar --infile',
            infile,
            '--img', imgfile,
            """--samps '{"samp1": ["int"]}'""",
            '--nterms 2 --meancol samp1_mean --target_rank genus',
            '--tabfile', tabfile,
        ])
        test_status = subprocess.call(cmd, shell=True)
        self.assertEqual(test_status, 0)
        nline = subprocess.run(['wc', '-l', tabfile], stdout=subprocess.PIPE)
        self.assertEqual(b'3', nline.stdout.strip().split()[0])
        os.remove(tabfile)

    def testFuncBar(self):
        infile = testfile('eggnog_out.tab')
        imgfile = testfile('test_eggnog_viz.png')
        samps = testfile('rudney_samples.tab')
        tabfile = testfile("eggnog_viz_file.tab")
        cmd = ' '.join([
            'python3 metaquantome/cli.py viz -m f --plottype bar '
            '--infile', infile,
            '--img', imgfile,
            '--samps', samps,
            '--nterms 20 --meancol NS_mean --target_onto bp',
            '--tabfile', tabfile
        ])
        test_status = subprocess.call(cmd, shell=True)
        self.assertEqual(test_status, 0)

    def testHeatmapViz(self):
        infile = testfile('ec_ttest_tested.tab')
        imgfile = testfile('cli_heatmap_viz.png')
        cmd = ' '.join([
            'python3 metaquantome/cli.py viz -m f --ontology ec --plottype heatmap',
            '--infile', infile,
            '--img', imgfile,
            "--samps '", TTEST_SINFO, "'",
            '--filter_to_sig', 'corrected_p_s1_over_s2',
            '--fc_corr_p', ,
            '--alpha 0.5'
        ])
        test_status = subprocess.call(cmd, shell=True)
        self.assertEqual(test_status, 0)
        cmd2 = ' '.join([
            'python3 metaquantome/cli.py viz -m f --ontology ec --plottype heatmap',
            '--infile', infile,
            '--img', imgfile,
            "--samps '", TTEST_SINFO, "'"
        ])
        test_status2 = subprocess.call(cmd2, shell=True)
        self.assertEqual(test_status2, 0)

    def testPCABig(self):
        infile = testfile('tax_filt_out.tab')
        sampfile = testfile('rudney_samples.tab')
        imgfile = testfile('cli_pca_viz.png')
        cmd = ' '.join([
            'python3 metaquantome/cli.py viz -m t --plottype pca',
            '--infile', infile,
            '--img', imgfile,
            "--samps", sampfile,
            '--calculate_sep'
        ])
        test_status = subprocess.call(cmd, shell=True)
        self.assertEqual(test_status, 0)


if __name__ == '__main__':
    unittest.main()
