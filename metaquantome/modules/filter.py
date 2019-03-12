import metaquantome.util.expand_io as expand_io
import metaquantome.util.stat_io as stat_io
from metaquantome.classes.SampleGroups import SampleGroups


def run_filter(expanded_file, sinfo, ontology, mode,
               qthreshold, min_child_non_leaf, min_child_nsamp, min_peptides,
               min_pep_nsamp, outfile=None):
    """
    Filter expanded dataframe to rows that satisfy filtering conditions.

    :param expanded_file: path to expanded file
    :param sinfo: Path to experimental design file
    :param ontology: relevant for f and ft modes. Either 'go', 'ec', or 'cog'
    :param mode: either 'f', 't', or 'ft'
    :param qthreshold: minimum number of quantitations per grp
    :param min_child_non_leaf: minimum number of children for terms that are not leaves
    :param min_child_nsamp: minimum number of samples with sample children greater than min_child_non_leaf
    :param min_peptides: minimum number of peptides for each term
    :param min_pep_nsamp: minimum number of samples where the number of peptides has to be larger than
    min_peptides
    :param outfile: File to write to
    :return: results dataframe
    """
    # create sample groups object
    samp_grps = SampleGroups(sinfo)

    # read in df
    df = stat_io.read_expanded_table(expanded_file, samp_grps)

    # filter to a minimum number of observed intensities per group
    samp_loc = samp_grps.grp_names.copy()

    # get the columns for the first sample group
    grp0 = samp_loc.pop(0)
    keep = get_rows_to_keep(mode, df, grp0, samp_grps, qthreshold=qthreshold, min_child_non_leaf=min_child_non_leaf,
                            min_child_nsamp=min_child_nsamp, min_peptides=min_peptides, min_pep_nsamp=min_pep_nsamp)

    # now, iterate over the rest
    while len(samp_loc) > 0:
        grp_i = samp_loc.pop(0)
        ith_keep = get_rows_to_keep(mode, df, grp_i, samp_grps, qthreshold=qthreshold,
                                    min_child_non_leaf=min_child_non_leaf, min_child_nsamp=min_child_nsamp,
                                    min_peptides=min_peptides, min_pep_nsamp=min_pep_nsamp)
        keep &= ith_keep
    filtered_df = df.loc[keep].copy()

    # write out
    if outfile:
        cols = expand_io.define_outfile_cols_expand(samp_grps, ontology, mode)
        expand_io.write_out_general(filtered_df, outfile, cols)
    return filtered_df


def get_rows_to_keep(mode, df, grp, samp_grps, qthreshold, min_child_non_leaf, min_child_nsamp, min_peptides,
                     min_pep_nsamp):
    """
    Use checking to find the rows (taxonomic or functional terms) that satisfy all of the filtering conditions for
    the specified group
    :param mode: either 'f', 't', or 'ft'
    :param df: data frame of functional and taxonomic terms. missing values are represented as 0.
    :param grp: grp to check conditions for
    :param samp_grps: SampleGroups() object
    :param qthreshold: minimum number of quantitations per grp
    :param min_child_non_leaf: minimum number of children for terms that are not leaves
    :param min_child_nsamp: minimum number of samples with sample children greater than min_child_non_leaf
    :param min_peptides: minimum number of peptides for each term
    :param min_pep_nsamp: minimum number of samples where the number of peptides has to be larger than
    min_peptides
    :return: boolean Series with rows to keep as True
    """
    # intensity
    intcols = samp_grps.sample_names[grp]
    keep_int = (df[intcols] > 0).apply(sum, axis=1) >= qthreshold

    # peptides
    peptide_cols = samp_grps.n_peptide_names_dict[grp]
    peptide_keep_series = (df[peptide_cols] > min_peptides)
    if min_pep_nsamp == "all":
        keep_peptide = peptide_keep_series.all(axis=1)
    else:
        keep_peptide = peptide_keep_series.apply(sum, axis=1) >= int(min_pep_nsamp)

    if mode != 'ft':
        # child non leaf
        child_cols = samp_grps.samp_children_names_dict[grp]
        child_keep_series = (df[child_cols] >= min_child_non_leaf) | (df[child_cols] == 0)
        if min_child_nsamp == "all":
            keep_child = child_keep_series.all(axis=1)
        else:
            keep_child = child_keep_series.apply(sum, axis=1) >= int(min_child_nsamp)
        all_keep = keep_int & keep_child & keep_peptide
    else:
        all_keep = keep_int & keep_peptide
    return all_keep
