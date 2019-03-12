class AnnotationNode:
    """
    A single taxon or functional term and its associated intensity.
    Also provides members for npeptide and n_sample_children
    """
    def __init__(self, id, intensity):
        """
        create AnnotationNode object

        :param id: unique id for the term
        :param intensity: a single number indicating the peptide intensity
        """
        self.id = id
        self.intensity = intensity
        self.npeptide = 1

        # the next two attributes are updated later in AnnotationHierarchy
        self.sample_children = None
        self.n_sample_children = None

    def add_peptide(self, intensity):
        """
        add peptide (aka term observation) to an existing
        node. Intensity is incremented by the new intensity,
        and npeptide is incremented by 1.

        :param intensity: observed peptide intensity associated with term
        :return: None
        """
        self.intensity += intensity
        if intensity > 0:
            self.npeptide += 1
