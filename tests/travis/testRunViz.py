import unittest
import os

from metaquantome.util.testutils import testfile, TTEST_SINFO
from metaquantome.analysis.run_viz import run_viz


class TestRunViz(unittest.TestCase):
    img=testfile('test.png')

    def testBasicTaxBar(self):
        infile = testfile('taxonomy_write_simple.tab')
        status = run_viz('bar', self.img, infile, 'tax',
                         nterms='2', meancol='samp1_mean',
                         target_rank="genus")
        self.assertEqual(status, 0)

    def testVolcano(self):
        infile = testfile('cli_mult_test_out.tab')
        status = run_viz('volcano', self.img, infile,
                         textannot="id",
                         fc_name="log2fc_s1_over_s2")
        self.assertEqual(status, 0)

    def testGOVolcano(self):
        infile = testfile('go_tested.tab')
        status = run_viz('volcano', self.img, infile,
                         textannot="id",
                         fc_name="log2fc_s1_over_s2",
                         gosplit=True)
        self.assertEqual(status, 0)

    def testHeatmap(self):
        infile = testfile('go_expanded_ttest.tab')
        status = run_viz('heatmap', self.img, infile,
                         sinfo=TTEST_SINFO)
        self.assertEqual(status, 0)

    def tearDown(self):
        if os.path.exists(self.img):
            os.remove(self.img)
