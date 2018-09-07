import unittest
from metaquant.AnnotationNode import AnnotationNode


class TestAnnotationNode(unittest.TestCase):
    def testInit(self):
        # simple creation
        pep1_intensity = [150, 100, 0]
        id = 'GO:0008150'
        an = AnnotationNode(id, pep1_intensity)
        self.assertEqual(an.intensity, pep1_intensity)
        self.assertEqual(an.id, id)
        self.assertEqual(an.npeptide, 1)

    def testAddPeptide(self):
        pep1_intensity = [150, 100, 0]
        id = 'GO:0008150'
        an = AnnotationNode(id, pep1_intensity)
        # add peptide evidence
        pep2_intensity = [100, 13, 32]
        expected = [x + y for x, y in zip(pep1_intensity, pep2_intensity)]
        an.add_peptide(pep2_intensity)
        self.assertEqual(an.intensity, expected)
        self.assertEqual(an.npeptide, 2)


