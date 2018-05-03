import unittest

import src.function as fq

class TestFunction(unittest.TestCase):
    def testSingleInt(self):
        go_df = fq.functional_analysis('go', 'data/test/function_single_int.tab', 'int', test=False)
        print(go_df)
        self.assertEqual(go_df.loc['GO:0008150']['int'], 1)
        self.assertEqual(go_df.loc["GO:0022610"]['int'], 2/3)
        self.assertEqual(go_df.loc["GO:0008152"]['int'], 1 / 3)

    def testMultipleInt(self):
        go_df = fq.functional_analysis('go', 'data/test/function_multiple_int.tab',
                                       grp1_intcols=['int1', 'int2', 'int3'], test=False)
        self.assertEqual(go_df.loc['GO:0008150']['int1'], 1)
        self.assertEqual(go_df.loc['GO:0008152']['int1'], 1/3)
        self.assertEqual(go_df.loc['GO:0022610']['int2'], 3/5)
        self.assertEqual(go_df.loc['GO:0008152']['int3'], 0.7)

    def testRedundant(self):
        """
        test that if a parent and child term are both present, the parent doesn't get double
        """
        godf = fq.functional_analysis('go',
                                      'data/test/function_multiple_go.tab',
                                      grp1_intcols=['int1', 'int2', 'int3'],
                                      test=False)
        self.assertEqual(godf.loc['GO:0008150']['int1'], 1)
        self.assertEqual(godf.loc['GO:0008152']['int1'], 1/3)
        self.assertEqual(godf.loc['GO:0022610']['int2'], 3/5)
        self.assertEqual(godf.loc['GO:0008152']['int3'], 0.7)

    def testEggnogOutput(self):
        go_df = fq.functional_analysis('go',
                                       'data/test/function_eggnog_gos.tabular',
                                       grp1_intcols=['int737WS', 'int737NS', 'int852WS',
                                                     'int852NS', 'int867WS', 'int867NS'],
                                       test=False)
        self.assertEqual(go_df.loc['GO:0008150']['int852WS'], 1)

        go_df = fq.functional_analysis('go', 'data/test/function_eggnog_gos.tabular',
                          grp1_intcols=['int737NS', 'int852NS', 'int867NS'],
                          grp2_intcols=['int737WS', 'int852WS', 'int867WS'],
                          test=True, slim_down=False, paired=True)

    def testDA(self):
        go_df = fq.functional_analysis('go', 'data/test/function_multiple_int_ttests.tab',
                                      grp1_intcols=['int1', 'int2', 'int3'],
                                      grp2_intcols=['int4', 'int5', 'int6'],
                                      test=True)
        # make sure all are less than 0.05
        self.assertTrue(go_df['corrected_p'].le(0.05).all())

    def testSlimDown(self):
        go_df = fq.functional_analysis('go', 'data/test/function_eggnog_gos.tabular',
                          grp1_intcols=['int737NS', 'int852NS', 'int867NS'],
                          grp2_intcols=['int737WS', 'int852WS', 'int867WS'],
                          test=True, slim_down=True, paired=True)

if __name__ == '__main__':
    unittest.main()