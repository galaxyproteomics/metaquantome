from metaquant.AnnotationHierarchy import AnnotationHierarchy
import pandas as pd


class SampleAnnotations:
    """
    Builds annotation hierarchy for each sample
    """
    def __init__(self, db):
        self.db = db
        self.hierarchies = set()
        self.df = None

    def add_sample(self, sample_set, df, annot_colname, int_colname):
        """
        Creates an AnnotationHierarchy for a given sample.
        :param sample_set: set of all observed annotations in that sample
        """
        samp = AnnotationHierarchy(self.db, sample_set, int_colname)
        # add node for each row in df
        samp.add_nodes_from_df(df, annot_colname, int_colname)
        return samp

    def add_samples_from_df(self, df, annot_colname, samp_grps, min_peptides, min_children_non_leaf):
        """
        Adds a sample for each intensity row in the dataframe
        :param df: Full dataframe
        :param annot_colname: Annotation column name
        :param samp_grps: SampleGroups object
        :return: Nothing
        """
        all_intcols = samp_grps.all_intcols
        hierarchies = set()
        for samp in all_intcols:
            filt = df.loc[df[samp] != 0]
            sample_set = set(df[annot_colname])
            hier = self.add_sample(sample_set, filt, annot_colname, samp)
            hier.get_informative_nodes(min_peptides, min_children_non_leaf)
            hierarchies.update({hier})
        self.hierarchies = hierarchies

    def to_dataframe(self):
        # convert each separate hierarchy to df
        # import pdb; pdb.set_trace()
        n_hier = len(self.hierarchies)
        loc_hier = self.hierarchies.copy()

        # remove any
        hierarchy_dfs = list()
        for i in range(n_hier):
            h = loc_hier.pop()
            if len(h.informative_nodes) > 0:
                dh = h.to_dataframe()
                hierarchy_dfs.append(dh)

        full_df = pd.concat(hierarchy_dfs, axis=1, sort=True)
        # stats expects 0's, not NaNs
        full_df.fillna(0)
        return full_df





