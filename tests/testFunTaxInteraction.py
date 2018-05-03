import unittest
import src.function_taxonomy_interaction as rp


class TestFunctionTaxInteraction(unittest.TestCase):
    def testReadAndDE(self):
        tf = rp.function_taxonomy_interaction_analysis('data/test/fun_tax_interaction.txt',
                                                       'cog', 'genus',
                                                       ['int1', 'int2', 'int3'],
                                                       ['int4', 'int5', 'int6'])
        cog1 = tf['function_taxon'][0]
        self.assertEqual(cog1, 'J-Streptococcus')
        self.assertTrue(all(tf['p.fdr'].le(0.05)))


if __name__=='__main__':
    unittest.main()