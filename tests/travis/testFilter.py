import unittest

from metaquantome.util.testutils import testfile, TTEST_SINFO
from metaquantome.util import stat_io
from metaquantome.util.utils import TEST_DIR
from metaquantome.classes.SampleGroups import SampleGroups
from metaquantome.modules.filter import run_filter
from metaquantome.modules.expand import expand


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

        expanded = expand('t', TTEST_SINFO, int_file=intfile, pep_colname_int='peptide', pep_colname_func='peptide',
                          pep_colname_tax='peptide', data_dir=TEST_DIR, outfile=expandfile, tax_file=taxfile,
                          tax_colname='lca')
        exp_ids = set(expanded['id'])

        # no filtering
        nofilt = run_filter(expandfile, TTEST_SINFO, ontology=None, mode="t", qthreshold=0, min_child_non_leaf=0,
                            min_child_nsamp=0, min_peptides=0, min_pep_nsamp=0)
        nofilt_ids = set(nofilt['id'])

        # make sure that ids are the same when no filtering is done
        self.assertSetEqual(nofilt_ids, exp_ids)

        # now, require 3 intensities per group. we shouldn't see 1496 or 1870884
        filt3 = run_filter(expandfile, TTEST_SINFO, ontology=None, mode="t", qthreshold=3, min_child_non_leaf=0,
                           min_child_nsamp=0, min_peptides=0, min_pep_nsamp=0)
        filt3_ids = set(filt3['id'])
        self.assertNotIn(1496, filt3_ids)
        self.assertNotIn(1870884, filt3_ids)


if __name__ == '__main__':
    unittest.main()
