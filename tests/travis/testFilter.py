import unittest

from metaquantome.util.testutils import testfile, TTEST_SINFO
from metaquantome.util import stat_io
from metaquantome.util.constants import TAX_TEST_DIR
from metaquantome.SampleGroups import SampleGroups
from metaquantome.analysis.filter import run_filter
from metaquantome.analysis.expand import expand


class TestFilter(unittest.TestCase):
    def testRead(self):
        samp_grps = SampleGroups('{"samp1": "int"}')
        infile = testfile('taxonomy_write_simple.tab')
        df = stat_io.read_expanded_table(infile, samp_grps)
        self.assertIn('Helicobacter', df['taxon_name'].tolist())

    def testFilter(self):
        intfile = testfile('filt_int.tab')
        taxfile = testfile('multiple_tax.tab')
        expandfile = testfile('expand_out.tab')

        expanded = expand('tax', TTEST_SINFO,
                          int_file = intfile,
                          tax_file=taxfile,
                          tax_colname='lca',
                          outfile=expandfile,
                          data_dir=TAX_TEST_DIR)
        exp_ids = set(expanded['id'])

        # no filtering
        nofilt = run_filter(expandfile, TTEST_SINFO,
                            qthreshold=0,
                            min_peptides=0,
                            min_pep_nsamp=0,
                            min_child_non_leaf=0,
                            min_child_nsamp=0,
                            mode="tax",
                            ontology=None)
        nofilt_ids = set(nofilt['id'])

        # make sure that ids are the same when no filtering is done
        self.assertSetEqual(nofilt_ids, exp_ids)

        # now, require 3 intensities per group. we shouldn't see 1496 or 1870884
        filt3 = run_filter(expandfile, TTEST_SINFO,
                           qthreshold=3,
                           min_peptides=0,
                           min_pep_nsamp=0,
                           min_child_non_leaf=0,
                           min_child_nsamp=0,
                           mode="tax",
                           ontology=None)
        filt3_ids = set(filt3['id'])
        self.assertNotIn(1496, filt3_ids)
        self.assertNotIn(1870884, filt3_ids)
