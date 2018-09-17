import numpy as np


class AnnotationNode:
    def __init__(self, id, intensity):
        """
        :param id: unique id for the term
        :param intensity: a list with intensity for each sample.
        The order will be kept constant by referring to the SampleGroups() object when calling this.
        """
        self.id = id
        self.intensity = intensity
        self.npeptide = self._calc_npeptide(intensity)

        # the next four attributes are updated later in AnnotationHierarchy
        self.sample_children = None
        self.sample_descendants = None
        self.n_sample_children = None
        self.aggregated_intensity = None

    def add_peptide(self, intensity):
        new_intensity = self._add_two_lists(self.intensity, intensity)
        self.intensity = new_intensity
        new_peptide_evidence = self._calc_npeptide(new_intensity)
        self.npeptide = self._add_two_lists(self.npeptide, new_peptide_evidence)

    def _calc_npeptide(self, intensity):
        return [(x > 0)*1 for x in intensity]

    def _add_two_lists(self, list1, list2):
        return [x + y for x, y in zip(list1, list2)]




