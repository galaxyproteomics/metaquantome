import unittest
from metaquant.SampleGroups import SampleGroups
from metaquant.util.stats import filter_min_observed
import pandas as pd


class TestStats(unittest.TestCase):
    def testFilterMinObserved(self):
        df = pd.DataFrame({'a1': [0, 1, 0, 1],
                           'a2': [1, 1, 0, 0],
                           'b1': [1, 1, 0, 0],
                           'b2': [1, 1, 1, 1]},
                          index=['pep1', 'pep2', 'pep3', 'pep4'])
        samp_grps = SampleGroups('{"a": ["a1", "a2"], "b": ["b1", "b2"]}')
        # first, if threshold is 0, should return original df
        filt0 = filter_min_observed(df, 0, samp_grps)
        self.assertTrue(df.equals(filt0))

        # if threshold is 1, should return pep1, pep2, and pep4
        filt1 = filter_min_observed(df, 1, samp_grps)
        exp1 = df.loc[['pep1', 'pep2', 'pep4']]
        self.assertTrue(exp1.equals(filt1))

        # if threshold is 2, should only return pep2
        filt2 = filter_min_observed(df, 2, samp_grps)
        exp2 = df.loc[['pep2']]
        self.assertTrue(exp2.equals(filt2))


if __name__ == '__main__':
    unittest.main()
