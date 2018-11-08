import unittest
import pandas as pd
import numpy as np

import metaquant.SampleGroups as sg
import metaquant.analysis.separation as sep


class TestSeparation(unittest.TestCase):
    sinfo = '{"a": ["a1", "a2"],"b": ["b1", "b2"]}'
    test_sg = sg.SampleGroups(sinfo)
    ind = ['t1', 't2']
    a1 = [100, 0]
    a2 = [20, 20]
    b1 = [0, 20]
    b2 = [1, 50]
    a_mean = np.array([60, 10])
    b_mean = np.array([0.5, 35])
    test_df = pd.DataFrame({'a1': a1,
                            'a2': a2,
                            'b1': b1,
                            'b2': b2},
                           index=ind)

    def testCentroid(self):
        test_means = sep.centroid(self.test_df, self.test_sg)
        mean_a = pd.Series([60, 10], index=self.ind)
        self.assertTrue((test_means[0] == mean_a).all())

    def testTotalSSEWithin(self):
        centroids = sep.centroid(self.test_df, self.test_sg)
        sse = sep.total_sse_within_group(self.test_df, centroids, self.test_sg)
        # expected sse
        tot_a = np.sum(np.square(self.a_mean - [100, 0])) +\
                np.sum(np.square(self.a_mean - [20, 20])) +\
                np.sum(np.square(self.b_mean - [0, 20])) +\
                np.sum(np.square(self.b_mean - [1, 50]))
        self.assertEqual(sse, tot_a)

    def testAvgSSEBetween(self):
        centroids = sep.centroid(self.test_df, self.test_sg)
        calc_dist = sep.avg_dist_between_group(centroids)
        expected_dist = np.sum(np.square(self.a_mean - self.b_mean))
        self.assertEqual(expected_dist, calc_dist)
