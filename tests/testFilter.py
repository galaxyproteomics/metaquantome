import unittest

from tests.testutils import testfile, TTEST_SINFO
from metaquant.util import io
from metaquant.SampleGroups import SampleGroups
from metaquant.analysis.filter import filter
from metaquant.analysis.expand import expand


class TestFilter(unittest.TestCase):
    def testRead(self):
        samp_grps = SampleGroups('{"samp1": "int"}')
        infile = testfile('taxonomy_write_simple.tab')
        df = io.read_expanded_table(infile, samp_grps)
        self.assertIn('Helicobacter', df['taxon_name'].tolist())

    def testFilter(self):
        intfile = testfile('filt_int.tab')
        taxfile = testfile('multiple_tax.tab')
        expandfile = testfile('expand_out.tab')

        expanded = expand('tax', TTEST_SINFO,
                          int_file = intfile,
                          tax_file=taxfile,
                          tax_colname='lca',
                          outfile=expandfile)
        exp_ids = set(expanded['id'])

        # no filtering
        nofilt = filter(expandfile, TTEST_SINFO,
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
        filt3 = filter(expandfile, TTEST_SINFO,
                       qthreshold=3,
                       min_peptides=0,
                       min_pep_nsamp=0,
                       min_child_non_leaf=0,
                       min_child_nsamp=0,
                       mode="tax",
                       ontology=None)
        print(filt3)
        filt3_ids = set(filt3['id'])
        self.assertNotIn(1496, filt3_ids)
        self.assertNotIn(1870884, filt3_ids)
        # os.remove(expandfile)
