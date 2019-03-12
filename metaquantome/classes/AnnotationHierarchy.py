import pandas as pd

import metaquantome.classes.AnnotationNode as anode
from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
from metaquantome.databases.NCBITaxonomyDb import NCBITaxonomyDb


class AnnotationHierarchy:
    """
    Annotation Hierarchy takes in a dataframe and builds a hierarchy for a specific sample.
    After hierarchy is built from the dataframe, the hierarchy can be written to dataframe.
    """
    def __init__(self, db, sample_set, sample_name):
        """
        Create instance of Annotation hierarchy and set methods.
        `nodes` and `informative_nodes` are filled later,
        in `add_nodes_from_df()`.

        :param db: The relevant database.
        :param sample_set: Set of all direct annotations in specific sample.
        :param sample_name: Name of the sample to build hierarchy for
        """
        self.db = db
        self.sample_set = sample_set
        self.nodes = dict()
        self.sample_name = sample_name

    def add_nodes_from_df(self, df, annot_colname, int_colname):
        """
        Adds all terms from a dataframe.

        :param df: Annotation dataframe. The dataframe should have one term per row.
        :param annot_colname: Name of annotation column.
        :param int_colname: Name of column with quantification for `self.sample_name`
        :return: None
        """
        for index, row in df.iterrows():
            term = row[annot_colname]
            # ensure each taxID is an integer (gets converted to float a lot by pandas)
            if isinstance(self.db, NCBITaxonomyDb):
                term = int(term)
            intensity = row[int_colname]
            self._add_node_with_ancestors(term, intensity)

        # add sample children here
        self._define_sample_children()

    def _add_node_with_ancestors(self, term, intensity):
        """
        create AnnotationNode objects for term and all of
        terms ancestors

        :param term: observed term
        :param intensity: intensity associated with term observation (peptide)
        :return: None
        """
        # add node for term itself
        self._add_node(term, intensity)

        # add node for ancestors
        ancestors = self.db.get_ancestors(term)
        for anc in ancestors:
            # if using GO and slimming down, only add ancestors in slim
            if isinstance(self.db, GeneOntologyDb):
                if self.db.slim_down and anc not in self.db.slim_members:
                    continue
            # self.expanded_sample_set.update({anc})
            self._add_node(anc, intensity)

    def _add_node(self, term, intensity):
        """
        create AnnotationNode for a given term observation
        or add intensity to an existing Annotation Node

        :param term: observed term
        :param intensity: intensity (numeric)
        :return: None
        """
        if term not in self.nodes.keys():
            # create new node
            self.nodes[term] = anode.AnnotationNode(term, intensity)
        else:
            # update existing node
            self.nodes[term].add_peptide(intensity)

    def _define_sample_children(self):
        """
        After all Node creation has happened, this
        defines the "sample children", or the children
        that each node has within the sample associated with
        the AnnotationHierarchy

        :return: None
        """
        # expanded sample set is all nodes
        expanded_sample_set = self.nodes.keys()

        for term in self.nodes.keys():
            node = self.nodes[term]
            ref_children = self.db.get_children(term)
            node.sample_children = ref_children.intersection(expanded_sample_set)
            node.n_sample_children = len(node.sample_children)

    def to_dataframe(self):
        """
        Run after add_nodes_from_df, and creates a dataframe
        that contains the aggregated intensity, number of sample children,
        and the number of peptides for each term

        :return: pandas data.frame with a row for each term
        """
        nodes = self.nodes
        node_rows = [0]*len(nodes)
        index = 0
        for term, node in nodes.items():
            samp_children_name = self.sample_name + '_n_samp_children'
            n_peptide_name = self.sample_name + '_n_peptide'
            node_rows[index] = pd.DataFrame({self.sample_name: node.intensity,
                                             samp_children_name: node.n_sample_children,
                                             n_peptide_name: node.npeptide}, index=[term])
            index += 1
        df = pd.concat(node_rows)
        return df
