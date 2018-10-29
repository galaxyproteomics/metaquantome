from metaquant.util import io
from metaquant.analysis.functional_analysis import functional_analysis
from metaquant.analysis.taxonomy_analysis import taxonomy_analysis
from metaquant.analysis.function_taxonomy_interaction import function_taxonomy_analysis
from metaquant.SampleGroups import SampleGroups

def expand(args):
    mode = args.mode
    samp_grps = SampleGroups(args.samps)

    # todo: move threshold to expand
    # read and join files - depending on pep/nopep
    df = io.read_and_join_files(mode, args.pep_colname, samp_grps,
                                int_file=args.int_file,
                                func_file=args.func_file, func_colname=args.func_colname,
                                tax_colname=args.tax_colname, tax_file=args.tax_file)
    # run analysis based on modes
    if mode == 'fn':
        results = functional_analysis(df=df, func_colname=args.func_colname, samp_grps=args.samps,
                                      ontology=args.ontology, slim_down=args.slim_down, data_dir=args.data_dir,
                                      overwrite=args.overwrite, min_peptides=args.min_peptides,
                                      min_children_non_leaf=args.min_children_non_leaf)
    elif mode == 'tax':
        results = taxonomy_analysis(df=df, samp_grps=samp_grps, data_dir=args.data_dir, tax_colname=args.tax_colname,
                                    min_peptides=args.min_peptides, min_children_non_leaf=args.min_children_non_leaf,
                                    threshold=args.threshold)
    elif mode == 'taxfn':
        if args.ontology != 'cog':
            raise ValueError("Only cog is supported for ft interaction. " +
                             "Make sure you have a cog column and supply the column name to func_colname")
        results = function_taxonomy_analysis(df=df, cog_name=args.func_colname, lca_colname=args.tax_colname,
                                             samp_grps=samp_grps,
                                             test=args.test, threshold=args.threshold,
                                             paired=args.paired, parametric=args.parametric,
                                             data_dir=args.data_dir)
    else:
        raise ValueError("Invalid mode. Expected one of: %s" % ['fun', 'tax', 'taxfn'])
    # filter min observed here

    # set up written output
    if args.outfile:
        cols = io.define_outfile_cols_expand(samp_grps, args.ontology, mode)
        results.to_csv(args.outfile, columns=cols, sep="\t", header=True, index=False, na_rep="NA")
    # whether writing out or not, return result data frame
    return results
