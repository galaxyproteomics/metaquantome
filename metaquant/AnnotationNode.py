# todo test both of these


class AnnotationNode:
    def __init__(self, id, intensity):
        """
        :param id: unique id for the term
        :param intensity: a list with intensity for each sample.
        The order will be kept constant by referring to the SampleGroups() object when calling this.
        """
        self.id = id
        self.intensity = intensity
        self.npeptide = 1

        # the next four attributes are updated later in AnnotationHierarchy
        self.sample_children = None
        self.sample_descendants = None
        self.n_sample_children = None
        self.aggregated_intensity = None

    def add_peptide(self, intensity):
        new_intensity = [x + y for x, y in zip(self.intensity, intensity)]
        self.intensity = new_intensity
        self.npeptide += 1


