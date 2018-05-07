import unittest
import metaquant
from definitions import DATA_DIR
import os

class TestFunctionTaxInteraction(unittest.TestCase):
    def testReadAndDE(self):
        datafile=os.path.join(DATA_DIR, 'test', 'fun_tax_interaction.txt')
        tf = metaquant.metaquant('taxfn',
                                 file=datafile,
                                 cog_colname='cog',
                                 tax_rank='genus',
                                 pep_colname="peptide",
                                 sample1_colnames=['int1', 'int2', 'int3'],
                                 sample2_colnames=['int4', 'int5', 'int6'])
        cog1 = tf['id'][0]
        self.assertEqual(cog1, 'J-Streptococcus')
        self.assertTrue(all(tf['corrected_p'].le(0.05)))


if __name__=='__main__':
    unittest.main()