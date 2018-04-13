import unittest
import diff_abund_tf.run_peca as rp


class TestReadIn(unittest.TestCase):
    def testRead(self):
        tf = rp.TaxaFuncPECA('test/test_read_tf.txt', 'cog', 'genus', ['int1', 'int2'], ['int3', 'int4'])
        cogs = tf.df.rx2('cog')
        self.assertEqual(cogs.levels[tf.df.rx2('cog')[0] - 1], 'M')


class TestPeca(unittest.TestCase):
    def testDe(self):
        tf = rp.TaxaFuncPECA('test/test_read_tf.txt', 'cog', 'genus',
                             ['int1', 'int2', 'int3'],
                             ['int4', 'int5', 'int6'])
        peca = tf.run_peca()
        self.assertTrue(all(peca['p.fdr'].le(0.05)))
