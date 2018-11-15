import pandas as pd

import metaquantome.AnnotationNode as anode
from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
from metaquantome.databases.NCBITaxonomyDb import NCBITaxonomyDb


class AnnotationHierarchy:
    """
    Annotation Hierarchy takes in a dataframe and builds a hierarchy for a specific sample.
    After hierarchy is build from the dataframe, the hierarchy can be written to dataframe.
    """
    def __init__(self, db, sample_set, sample_name):
        """
        Create instance of Annotation hierarchy and set methods.
        Some are filled later, in `add_nodes_from_df()`. These are
        `nodes` and `informative_nodes`
        :param db: The relevant database.
        :param sample_set: Set of all direct annotations in specific sample.
        :param sample_name: Name of the sample to build hierarchy for
        """
        self.db = db
        self.sample_set = sample_set
        self.expanded_sample_set = sample_set.copy()  # add more terms later
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
            if isinstance(self.db, NCBITaxonomyDb):  # todo: move this to IO
                term = int(term)
            intensity = row[int_colname]
            self._add_node(term, intensity)
        # add sample children here
        self._define_sample_children()

    def to_dataframe(self):
        # todo: doc
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

    def _add_node(self, term, intensity):
        # todo: doc
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
            self._add_node(par, intensity)

    def _define_sample_children(self):
        # todo: doc
        for term in self.nodes.keys():
            node = self.nodes[term]
            ref_children = self.db.get_children(term)
            node.sample_children = ref_children.intersection(self.expanded_sample_set)
            node.n_sample_children = len(node.sample_children)
