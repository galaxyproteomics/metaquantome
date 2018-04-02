import unittest
import func.goquant as gq
import pandas as pd


class TestSingleInt(unittest.TestCase):
    def testBioProc(self):
        # both are immediate children of biological_process (GO:0008150)
        go1 = "GO:0008152"
        go2 = "GO:0022610"
        pep = ["A", "B"]
        inten = [10, 20]
        df = pd.DataFrame({"go": [go1, go2], "int": inten, "pep": pep})
        godf = gq.GoQuant(df, 'int')
        print(godf.go_int)
        self.assertEqual(godf.go_int.loc['GO:0008150']['intensity'], 30)
        self.assertEqual(godf.go_int.loc[go2]['OSRA'], 2/3)
        self.assertEqual(godf.go_int.loc[go1]['OSRA'], 1 / 3)


class TestMultInt(unittest.TestCase):
    def testThree(self):
        go1 = "GO:0008152"
        go2 = "GO:0022610"
        pep = ["A", "B"]
        int1 = [10, 20]
        int2 = [20, 30]
        int3 = [70, 30]
        df = pd.DataFrame({"go": [go1, go2], "int1": int1,
                           "int2": int2, "int3": int3, "pep": pep})
        godf = gq.GoQuant(df, ['int1', 'int2', 'int3'])
        self.assertEqual(godf.go_int.loc['GO:0008150']['OSRA_int1'], 1)
        self.assertEqual(godf.go_int.loc['GO:0008152']['OSRA_int1'], 1/3)
        self.assertEqual(godf.go_int.loc['GO:0022610']['OSRA_int2'], 3/5)
        self.assertEqual(godf.go_int.loc['GO:0008152']['OSRA_int3'], 0.7)

if __name__ == '__main__':
    unittest.main()