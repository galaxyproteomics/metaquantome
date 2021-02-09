import unittest
import os
import pandas as pd

from metaquantome.util.testutils import testfile, TTEST_SINFO
from metaquantome.modules.run_viz import run_viz


class TestRunViz(unittest.TestCase):
    img = testfile('test.png')

    def testBasicTaxBar(self):
        infile = testfile('taxonomy_write_simple.tab')
        tabfile = testfile('taxonomy_plot_out.tab')
        run_viz('bar', self.img, infile, mode='t',
                nterms='2', meancol='samp1_mean',
                target_rank="genus", barcol="6",
                tabfile=tabfile)
        self.assertTrue(os.path.exists(tabfile))
        os.remove(tabfile)
        run_viz('bar', self.img, infile, mode='t',
                nterms='2', meancol='samp1_mean',
                target_rank="genus", barcol="6")

    def testFuncBar(self):
        infile = testfile('eggnog_out.tab')
        tabfile = testfile('viz_f_filt.tab')
        # bp
        run_viz('bar', self.img, infile, mode='f',
                nterms='2', meancol='NS_mean',
                target_onto='bp', barcol='1', tabfile=tabfile)
        df = pd.read_csv(tabfile, sep='\t')
        namespace = df.namespace.unique()
        self.assertEqual(namespace, ['biological_process'])

        # cc
        run_viz('bar', self.img, infile, mode='f',
                nterms='2', meancol='NS_mean',
                target_onto='cc', barcol='1', tabfile=tabfile)
        df = pd.read_csv(tabfile, sep='\t')
        namespace = df.namespace.unique()
        self.assertEqual(namespace, ['cellular_component'])

        # mf
        run_viz('bar', self.img, infile, mode='f',
                nterms='2', meancol='NS_mean',
                target_onto='mf', barcol='1', tabfile=tabfile)
        df = pd.read_csv(tabfile, sep='\t')
        namespace = df.namespace.unique()
        self.assertEqual(namespace, ['molecular_function'])

    def testVolcano(self):
        infile = testfile('cli_mult_test_out.tab')
        run_viz('volcano', self.img, infile,
                textannot="id",
                fc_corr_p="corrected_p_s1_over_s2",
                fc_name="log2fc_s1_over_s2")

        # test tabfile
        tabfile = testfile('taxonomy_plot_out.tab')
        run_viz('volcano', self.img, infile,
                textannot="id",
                fc_name="log2fc_s1_over_s2",
                fc_corr_p="corrected_p_s1_over_s2",
                tabfile=tabfile)
        self.assertTrue(os.path.exists(tabfile))
        os.remove(tabfile)

    def testGOVolcano(self):
        infile = testfile('go_tested.tab')
        run_viz('volcano', self.img, infile,
                textannot="id",
                fc_name="log2fc_s1_over_s2",
                fc_corr_p="corrected_p_s1_over_s2",
                gosplit=True)

    def testHeatmap(self):
        infile = testfile('go_expanded_ttest.tab')
        run_viz('heatmap', self.img, infile,
                sinfo=TTEST_SINFO)

    def testPCA(self):
        infile = testfile('go_tested.tab')
        run_viz('pca', self.img, infile,
                sinfo=TTEST_SINFO,
                calculate_sep=False)

    def testPCABig(self):
        infile = testfile('tax_filt_out.tab')
        sampfile = testfile('rudney_samples.tab')
        run_viz('pca', self.img, infile,
                sinfo=sampfile,
                calculate_sep=True)

    def testFtDist(self):
        infile = testfile('ft_out.tab')
        run_viz('ft_dist', self.img, infile,
                meancol="s1_mean",
                whichway='t_dist',
                id="GO:0008150",
                target_rank="genus",
                nterms="all")

        # f dist
        run_viz('ft_dist', self.img, infile,
                meancol="s1_mean",
                whichway='f_dist',
                target_onto="bp",
                nterms="all", id=209)

        # test tabfile
        tabfile = testfile("tmp")
        run_viz('ft_dist', self.img, infile,
                meancol="s1_mean",
                whichway='t_dist',
                id="GO:0008150",
                target_rank="genus",
                nterms="all",
                tabfile=tabfile)
        self.assertTrue(os.path.exists(tabfile))
        os.remove(tabfile)

    # comment this out for test plot inspection
    def tearDown(self):
        if os.path.exists(self.img):
            os.remove(self.img)


if __name__=='__main__':
    unittest.main()
