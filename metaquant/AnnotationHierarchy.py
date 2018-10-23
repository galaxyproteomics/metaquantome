import metaquant.AnnotationNode as anode
import pandas as pd
from metaquant.databases.GeneOntologyDb import GeneOntologyDb
from metaquant.databases.NCBITaxonomyDb import NCBITaxonomyDb
from metaquant.util import utils

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
    def __init__(self, db, sample_set, sample_name):
        self.db = db
        self.sample_set = sample_set
        self.expanded_sample_set = sample_set  # add more terms later
        self.nodes = dict()
        self.informative_nodes = dict()
        self.sample_name = sample_name

    def add_nodes_from_df(self, df, annot_colname, int_colname):
        for index, row in df.iterrows():
            term = row[annot_colname]
            if isinstance(self.db, NCBITaxonomyDb):  # todo: move this to IO
                term = int(term)
            intensity = row[int_colname]
            self.add_node(term, intensity)

    def add_node(self, term, intensity):
        if term not in self.nodes.keys():
            # create new node
            self.nodes[term] = anode.AnnotationNode(term, intensity)
        else:
            # update existing node
            self.nodes[term].add_peptide(intensity)
        # do same for parents
        parents = self.db.get_parents(term)
        for par in parents:
            # if using GO and slimming down, only add parents in slim
            if isinstance(self.db, GeneOntologyDb):
                if self.db.slim_down and par not in self.db.slim_members:
                    continue
            self.expanded_sample_set.update({par})
            self.add_node(par, intensity)

    def _define_sample_children(self):
        for term in self.nodes.keys():
            node = self.nodes[term]
            ref_children = self.db.get_children(term)
            node.sample_children = ref_children.intersection(self.expanded_sample_set)
            node.n_sample_children = len(node.sample_children)

    def get_informative_nodes(self, min_peptides, min_children_non_leaf):
        """
        :param min_peptides: if node has fewer than min_peptides (i.e. < min_peptides), will be pruned
        :param min_children_non_leaf: if node has fewer (<) than min_children_non_leaf and more than 0 children,
        it will be pruned.
        :return:
        """
        # first, define sample children for each term
        self._define_sample_children()
        informative_nodes = dict()
        for term, node in self.nodes.items():
            n_children = node.n_sample_children
            n_peptides = node.npeptide
            if (n_children >= min_children_non_leaf or n_children == 0) and n_peptides >= min_peptides:
                # add node to informative node dict
                informative_nodes[term] = node
        # change self nodes to informative nodes
        self.informative_nodes = informative_nodes

    def to_dataframe(self):
        inf_nodes = self.informative_nodes
        node_rows = [0]*len(inf_nodes)
        index = 0
        for term, node in inf_nodes.items():
            node_rows[index] = pd.DataFrame({self.sample_name: node.intensity}, index=[term])
            index += 1
        df = pd.concat(node_rows)
        return df
