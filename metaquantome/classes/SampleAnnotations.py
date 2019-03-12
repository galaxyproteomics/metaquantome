import pandas as pd

from metaquantome.classes.AnnotationHierarchy import AnnotationHierarchy


class SampleAnnotations:
    """
    Builds annotation hierarchy for each sample by
    calling `AnnotationHierarchy`
    """
    def __init__(self, db):
        """
        Initialize SampleAnnotations object

        :param db: reference database - will be taxonomy, GO terms, etc.
        """
        self.db = db
        # this is filled later in add_samples_from_df
        self.hierarchies = set()

    def add_samples_from_df(self, df, annot_colname, samp_grps):
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
            filt = self.filter_to_samp_observed(df, samp)
            sample_set = self.create_sample_set(filt, annot_colname)
            # create hierarchy for this sample
            hier = AnnotationHierarchy(self.db, sample_set, samp)
            # add node for each row in df
            hier.add_nodes_from_df(filt, annot_colname, samp)
            # add to set
            hierarchies.update({hier})
        self.hierarchies = hierarchies

    @staticmethod
    def create_sample_set(df_filt, annot_colname):
        """
        Create the set of all terms in the specific sample

        :param df_filt: DataFrame with one term per row, filtered so that only terms present in the sample of interest
        are included in the annot_colname column
        :param annot_colname: Name of annotation column - either tax_colname or func_colname
        :return: Set of all terms in df_filt
        """
        return set(df_filt[annot_colname])

    @staticmethod
    def filter_to_samp_observed(df, samp):
        # todo: add test
        """
        Filter to observed values in specific sample.

        :param df: DataFrame with all intensities and annotations. Missing values should be 0.
        :param samp: The specific sample of interest
        :return: DataFrame that only has rows whether the quantified value for the sample of interest is
        nonzero.
        """
        filt = df.loc[df[samp] != 0]
        return filt

    def to_dataframe(self):
        """
        convert each hierarchy to dataframe and concatenate.
        fills NAs with 0s, for adding up
        :return: concatenated dataframe, missing values are 0
        """
        n_hier = len(self.hierarchies)
        loc_hier = self.hierarchies.copy()
        hierarchy_dfs = list()
        for i in range(n_hier):
            h = loc_hier.pop()
            dh = h.to_dataframe()
            hierarchy_dfs.append(dh)
        full_df = pd.concat(hierarchy_dfs, axis=1, sort=True)
        # stats expects 0's, not NaNs
        full_df.fillna(0)
        return full_df
