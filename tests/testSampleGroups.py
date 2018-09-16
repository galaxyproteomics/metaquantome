import unittest

import metaquant.SampleGroups
from metaquant import io
import os
from metaquant.utils import DATA_DIR


class TestSampleGroups(unittest.TestCase):

    def testReadSampInfoFile(self):
        # test that info file is correctly parsed
        sinfo = os.path.join(DATA_DIR, 'test', 'samp_info.tab')
        sample_groups = metaquant.SampleGroups.SampleGroups(sinfo)
        self.assertDictEqual(sample_groups.sample_names, {'A': ['A1', 'A2'], 'B': ['B1', 'B2']})

    def testReadSampInfoJson(self):
        sinfo = '{"A": ["A1", "A2"], "B": ["B1", "B2"]}'
        samp_grps = metaquant.SampleGroups.SampleGroups(sinfo)
        self.assertEqual(samp_grps.grp_names, ['A', 'B'])

        # note that the order of all_intcols is random (so we compare sets)
        self.assertEqual(set(samp_grps.all_intcols), {'A1', 'A2', 'B1', 'B2'})

if __name__ == '__main__':
    unittest.main()