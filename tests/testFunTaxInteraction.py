import unittest
import metaquant
from definitions import DATA_DIR
import os

class TestFunctionTaxInteraction(unittest.TestCase):
    def testReadAndDE(self):
        datafile=os.path.join(DATA_DIR, 'test', 'fun_tax_interaction.txt')
        outfile=os.path.join(DATA_DIR, 'test', 'fun_tax_interaction_out.txt')
        tf = metaquant.metaquant('taxfn',
                                 file=datafile,
                                 cog_colname='cog',
                                 lca_colname='genus',
                                 pep_colname="peptide",
                                 test=True,
                                 sample_names={'samp1': ['int1', 'int2', 'int3'],
                                               'samp2': ['int4', 'int5', 'int6']},
                                 outfile=outfile)
        cog1 = tf['id'][0]
        self.assertEqual(cog1, 'J-Streptococcus')
        self.assertTrue(all(tf['corrected_p'].le(0.05)))
        # check that ratios are in right direction: M-Streptococcus should be up and J-Strep should be down
        self.assertTrue(tf.loc['M-Streptococcus']['log2fc_samp1_over_samp2'] > 0)
        self.assertTrue(tf.loc['J-Streptococcus']['log2fc_samp1_over_samp2'] < 0)


if __name__=='__main__':
    unittest.main()