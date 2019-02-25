import unittest

from metaquantome.classes.AnnotationNode import AnnotationNode


class TestAnnotationNode(unittest.TestCase):
    def testInit(self):
        # simple creation
        pep1_intensity = 150
        id = 'GO:0008150'
        an = AnnotationNode(id, pep1_intensity)
        self.assertEqual(an.intensity, 150)
        self.assertEqual(an.id, id)
        self.assertEqual(an.npeptide, 1)

    def testAddPeptide(self):
        pep1_intensity = 100
        id = 'GO:0008150'
        an = AnnotationNode(id, pep1_intensity)
        # add peptide evidence
        pep2_intensity = 0
        an.add_peptide(pep2_intensity)
        self.assertEqual(an.intensity, pep1_intensity)
        self.assertEqual(an.npeptide, 1)


