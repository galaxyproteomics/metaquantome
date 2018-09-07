import metaquant.AnnotationNode as anode

# Annotation Hierarchy that takes in the dataframe and builds hierarchy
# the pruning method removes all nodes with a number of children less than N
# and a number of peptides less than M
# but the total number of peptides for a node depends on the sum of the peptides of the node's descendants
# so create a total_peptides class member that calculates the peptides of all descendants


class AnnotationHierarchy:
    """
    in init, define properties:
        - nodes: empty dict
    then, in update_node,
        1) if node doesn't already exist, create it.
        2) call add_peptide method on the relevant node
    methods:
        - prune(): filter all nodes based on evidence and informativeness
        - to_dataframe() for collapsing to dataframe
    """
    def __init__(self, db, sample_set):
        self.db = db
        self.sample_set = sample_set
        self.nodes = dict()
        self.informative_nodes = dict()

    def update_node(self, id, intensity):
        if id not in self.nodes.keys():
            # create new node
            new_node = anode.AnnotationNode(id, intensity)

            # children handling #
            # get node's children from reference database
            ref_children = self.db.get_children(id)
            # get node's children that are also present in sample_set, and assign to this node
            new_node.sample_children = ref_children.intersection(self.sample_set)

            # determine number of sample children
            new_node.n_sample_children = len(new_node.sample_children)

            # descendants handling #
            ref_descendants = self.db.get_descendants(id)
            new_node.sample_descendants = ref_descendants.intersection(self.sample_set)

            # insert into dict
            self.nodes[id] = new_node

        else:
            # update existing node
            self.nodes[id].add_peptide(intensity)

    def aggregate_nodes(self):
        # this sums up intensity for all nodes:
        # get sample descendants
        # sum intensities of all sample descendants
        # something like [sum(x) for x in zip(*intensities)]
        for id, node in self.nodes.items():
            descendant_ids = node.sample_descendants
            intensities = [node.intensity]
            for id in descendant_ids:
                descendant_node = self.nodes[id]
                intensities.append(descendant_node.intensity)
            # aggregate intensity
            agg_intensity = [sum(x) for x in zip(*intensities)]
            node.aggregated_intensity = agg_intensity  # make sure this updates the intensity in the node itself

    def get_informative_nodes(self, min_peptides, min_children_non_leaf):
        """
        :param min_peptides: if node has fewer than min_peptides (i.e. < min_peptides), will be pruned
        :param min_children_non_leaf: if node has fewer (<) than min_children_non_leaf and more than 0 children,
        it will be pruned.
        :return:
        """
        informative_nodes = dict()
        for id, node in self.nodes.items():
            n_children = node.n_sample_children
            n_peptides = node.npeptide
            if (n_children >= min_children_non_leaf or n_children == 0) and n_peptides >= min_peptides:
                # add node to informative node dict
                informative_nodes[id] = node
        return informative_nodes

    def to_dataframe(self):
        # todo
        pass