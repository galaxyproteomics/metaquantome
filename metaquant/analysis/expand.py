from metaquant.util import io
from metaquant.analysis.functional_analysis import functional_analysis
from metaquant.analysis.taxonomy_analysis import taxonomy_analysis
from metaquant.analysis.function_taxonomy_interaction import function_taxonomy_analysis
from metaquant.SampleGroups import SampleGroups


def expand(mode, samps, int_file, pep_colname='peptide', data_dir=None, overwrite=False, outfile=None, func_file=None,
           func_colname=None, ontology='go', slim_down=False, tax_file=None, tax_colname=None, min_peptides=0,
           min_children_non_leaf=0, threshold=0):

    samp_grps = SampleGroups(samps)

    # todo: move threshold to expand
    # read and join files - depending on pep/nopep
    df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                int_file=int_file,
                                func_file=func_file, func_colname=func_colname,
                                tax_colname=tax_colname, tax_file=tax_file)
    # run analysis based on modes
    if mode == 'fn':
        results = functional_analysis(df=df, func_colname=func_colname, samp_grps=samp_grps,
                                      ontology=ontology, slim_down=slim_down, data_dir=data_dir,
                                      overwrite=overwrite, min_peptides=min_peptides,
                                      min_children_non_leaf=min_children_non_leaf)
    elif mode == 'tax':
        results = taxonomy_analysis(df=df, samp_grps=samp_grps, data_dir=data_dir, tax_colname=tax_colname,
                                    min_peptides=min_peptides, min_children_non_leaf=min_children_non_leaf,
                                    threshold=threshold)
    elif mode == 'taxfn':
        # if ontology != 'cog':
        #     raise ValueError("Only cog is supported for ft interaction. " +
        #                      "Make sure you have a cog column and supply the column name to func_colname")
        results = function_taxonomy_analysis(df=df, cog_name=func_colname, lca_colname=tax_colname, samp_grps=samp_grps,
                                             threshold=threshold, data_dir=data_dir)
    else:
        raise ValueError("Invalid mode. Expected one of: %s" % ['fun', 'tax', 'taxfn'])
    # filter min observed here

    # set up written output
    if outfile:
        cols = io.define_outfile_cols_expand(samp_grps, ontology, mode)
        results.to_csv(outfile, columns=cols, sep="\t", header=True, index=False, na_rep="NA")
    # whether writing out or not, return result data frame
    return results
