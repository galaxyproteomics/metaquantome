import metaquantome.util.expand_io as expand_io
import metaquantome.util.stat_io as stat_io
from metaquantome.SampleGroups import SampleGroups


def run_filter(file, sinfo, ontology, mode,
               qthreshold,
               min_child_non_leaf, min_child_nsamp,
               min_peptides, min_pep_nsamp, outfile=None):
    # todo: doc
    # create sample groups object
    samp_grps = SampleGroups(sinfo)

    # read in df
    df = stat_io.read_expanded_table(file, samp_grps)

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

    if mode != 'taxfn':
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
