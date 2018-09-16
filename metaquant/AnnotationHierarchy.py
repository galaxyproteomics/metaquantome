import metaquant.AnnotationNode as anode
import pandas as pd
from metaquant.GeneOntologyDb import GeneOntologyDb
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
        self.expanded_sample_set = sample_set  # add more terms later
        self.nodes = dict()
        self.informative_nodes = dict()

    def add_node(self, id, intensity):
        if id not in self.nodes.keys():
            # create new node
            self.nodes[id] = anode.AnnotationNode(id, intensity)
        else:
            # update existing node
            self.nodes[id].add_peptide(intensity)
        # do same for parents #
        parents = self.db.get_parents(id)
        for par in parents:
            # if using GO and slimming down, only add parents in slim
            if isinstance(self.db, GeneOntologyDb):
                if self.db.slim_down and par not in self.db.slim_members:
                    continue
            self.expanded_sample_set.update({par})
            self.add_node(par, intensity)

    def _define_sample_children(self):
        for id in self.nodes.keys():
            node = self.nodes[id]
            ref_children = self.db.get_children(id)
            node.sample_children = ref_children.intersection(self.expanded_sample_set)
            node.n_sample_children = len(node.sample_children)

    # def aggregate_nodes(self):
    #     # this sums up intensity for all nodes:
    #     # get sample descendants
    #     # sum intensities of all sample descendants
    #     # something like [sum(x) for x in zip(*intensities)]
    #     for id, node in self.nodes.items():
    #         descendant_ids = node.sample_descendants
    #         intensities = [node.intensity]
    #         for id in descendant_ids:
    #             descendant_node = self.nodes[id]
    #             intensities.append(descendant_node.intensity)
    #         # aggregate intensity
    #         agg_intensity = [sum(x) for x in zip(*intensities)]
    #         node.aggregated_intensity = agg_intensity  # make sure this updates the intensity in the node itself

    def get_informative_nodes(self, min_peptides, min_children_non_leaf):
        """
        :param min_peptides: if node has fewer than min_peptides (i.e. < min_peptides), will be pruned
        :param min_children_non_leaf: if node has fewer (<) than min_children_non_leaf and more than 0 children,
        it will be pruned.
        :return:
        """
        # first, define sample children for each id
        self._define_sample_children()
        informative_nodes = dict()
        for id, node in self.nodes.items():
            n_children = node.n_sample_children
            n_peptides = node.npeptide
            if (n_children >= min_children_non_leaf or n_children == 0) and n_peptides >= min_peptides:
                # add node to informative node dict
                informative_nodes[id] = node
        # change self nodes to informative nodes
        self.nodes = informative_nodes

    def to_dataframe(self, intensity_names):
        node_rows = [0]*len(self.nodes)
        index = 0
        for uid, node in self.nodes.items():
            intensity = node.intensity
            # put aggregated intensities into dict
            ints = {name: value for name, value in zip(intensity_names, intensity)}
            node_rows[index] = pd.DataFrame(ints, index=[uid])
            index += 1
        df = pd.concat(node_rows)
        df['id'] = df.index
        return df
