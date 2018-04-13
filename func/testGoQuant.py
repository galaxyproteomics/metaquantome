import unittest
import func.goquant as gq
import pandas as pd


class TestSingleInt(unittest.TestCase):
    def testBioProc(self):
        godf = gq.GoQuant('int', 'go', 'test/sing_int.tab')
        self.assertEqual(godf.go_int.loc['GO:0008150']['int'], 30)
        self.assertEqual(godf.go_int.loc["GO:0022610"]['OSRA_int'], 2/3)
        self.assertEqual(godf.go_int.loc["GO:0008152"]['OSRA_int'], 1 / 3)


class TestMultInt(unittest.TestCase):
    def testThree(self):
        godf = gq.GoQuant(['int1', 'int2', 'int3'], 'go', 'test/mult_int.tab')
        self.assertEqual(godf.go_int.loc['GO:0008150']['OSRA_int1'], 1)
        self.assertEqual(godf.go_int.loc['GO:0008152']['OSRA_int1'], 1/3)
        self.assertEqual(godf.go_int.loc['GO:0022610']['OSRA_int2'], 3/5)
        self.assertEqual(godf.go_int.loc['GO:0008152']['OSRA_int3'], 0.7)

    def testRedundant(self):
        """
        test that if a parent and child term are both present, the parent doesn't get double
        """
        godf = gq.GoQuant(['int1', 'int2', 'int3'], 'go', 'test/mult_go.tab')
        self.assertEqual(godf.go_int.loc['GO:0008150']['OSRA_int1'], 1)
        self.assertEqual(godf.go_int.loc['GO:0008152']['OSRA_int1'], 1/3)
        self.assertEqual(godf.go_int.loc['GO:0022610']['OSRA_int2'], 3/5)
        self.assertEqual(godf.go_int.loc['GO:0008152']['OSRA_int3'], 0.7)


class TestRealData(unittest.TestCase):
    def testEggnog(self):
        godf = gq.GoQuant(['int737WS', 'int737NS', 'int852WS', 'int852NS', 'int867WS', 'int867NS'], 'go', 'test/eggnog_gos.tabular')
        self.assertEqual(godf.go_int.loc['GO:0008150']['OSRA_int852WS'], 1)

if __name__ == '__main__':
    unittest.main()