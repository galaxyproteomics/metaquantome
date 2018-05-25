import unittest
import numpy as np
import pandas as pd
from src import common


class TestCommon(unittest.TestCase):
    def testDefineIntensityColumns(self):
        # single string
        sing = {'y': ['y1']}
        samps = common.SampleGroups(sing)
        self.assertEqual(samps.all_intcols, ['y1'])
        self.assertEqual(samps.dict_numeric_cols, {x : np.float64 for x in samps.all_intcols})

        sample_names = {'s': ['s1', 's2', 's3'],
                        'k': ['k1', 'k2', 'k3'],
                        'x': ['x1', 'x2', 'x3'],
                        'y': ['y1']}

        # the order of dictionaries is not defined, so compare sets
        flattened_samps = {'s1', 's2', 's3',
                           'k1', 'k2', 'k3',
                           'x1', 'x2', 'x3',
                           'y1'}

        samps = common.SampleGroups(sample_names)
        self.assertEqual(flattened_samps, set(samps.all_intcols))
        self.assertEqual(samps.dict_numeric_cols, {x: np.float64 for x in flattened_samps})

    def testCalcMeans(self):
        df = pd.DataFrame({'s1_1': [4, 4],
                           's1_2': [2, 2],
                           's2_1': [5, 10],
                           's2_2': [7, 16]})
        samps = common.SampleGroups({'s1': ['s1_1', 's1_2'], 's2': ['s2_1', 's2_2']})
        means = common.calc_means(df, samps)
        self.assertTrue(means['s1_mean'].equals(pd.Series({0: 3.0, 1:3.0}, name="s1_mean")))

    def testFoldChange(self):
        df = pd.DataFrame({'s1_1': [4, 4],
                           's1_2': [2, 2],
                           's2_1': [5, 10],
                           's2_2': [7, 16]})
        samps = common.SampleGroups({'s1': ['s1_1', 's1_2'], 's2': ['s2_1', 's2_2']})

        means = common.calc_means(df, samps)
        fc = common.fold_change(means, samps)
        self.assertTrue(fc['log2fc_s1_over_s2'].equals(pd.Series({0: np.log2(3/6), 1: np.log2(3/13)})))

if __name__=='__main__':
    unittest.main()