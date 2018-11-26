import unittest
import os

from metaquantome.util.testutils import testfile
from metaquantome.analysis.run_viz import run_viz


class TestRunViz(unittest.TestCase):
    def testBasicTax(self):
        infile = testfile('taxonomy_write_simple.tab')
        status = run_viz('bar', 'test.png', infile, 'tax',
                nterms='2', meancol='samp1_mean',
                target_rank="genus")
        self.assertEqual(status, 0)
        os.remove('test.png')
