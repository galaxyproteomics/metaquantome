import unittest
import os
import numpy as np

import metaquantome.classes.SampleGroups
from metaquantome.util import utils


class TestSampleGroups(unittest.TestCase):

    def testReadSampInfoFile(self):
        # test that info file is correctly parsed
        sinfo = os.path.join(utils.DATA_DIR, 'test', 'samp_info.tab')
        sample_groups = metaquantome.classes.SampleGroups.SampleGroups(sinfo)
        self.assertDictEqual(sample_groups.sample_names, {'A': ['A1', 'A2'], 'B': ['B1', 'B2']})

    def testReadSampInfoJson(self):
        sinfo = '{"A": ["A1", "A2"], "B": ["B1", "B2"]}'
        samp_grps = metaquantome.classes.SampleGroups.SampleGroups(sinfo)
        self.assertEqual(samp_grps.grp_names, ['A', 'B'])

        # note that the order of all_intcols is random (so we compare sets)
        self.assertEqual(samp_grps.all_intcols, ['A1', 'A2', 'B1', 'B2'])

    def testDefineIntensityColumns(self):
        # single string
        sing = '{"y": ["y1"]}'
        samps = metaquantome.classes.SampleGroups.SampleGroups(sing)
        self.assertEqual(samps.all_intcols, ['y1'])
        self.assertEqual(samps.dict_numeric_cols, {x: np.float64 for x in samps.all_intcols})

        sample_names = '{"s": ["s1", "s2", "s3"],"k": ["k1", "k2", "k3"],"x": ["x1", "x2", "x3"],"y": ["y1"]}'

        # the order of dictionaries is not defined, so compare sets
        flattened_samps = {'s1', 's2', 's3',
                           'k1', 'k2', 'k3',
                           'x1', 'x2', 'x3',
                           'y1'}

        samps = metaquantome.classes.SampleGroups.SampleGroups(sample_names)
        self.assertEqual(flattened_samps, set(samps.all_intcols))
        self.assertEqual(samps.dict_numeric_cols, {x: np.float64 for x in flattened_samps})


if __name__ == '__main__':
    unittest.main()
