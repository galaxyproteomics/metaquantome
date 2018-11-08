import sys
import argparse
import logging


from metaquantome.analysis.expand import expand
from metaquantome.analysis.filter import run_filter
from metaquantome.analysis.stat import stat


def cli():
    # initialize logger
    logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stderr)
    args = parse_args_cli()
    if args.command == "expand":
        expand(mode=args.mode, samps=args.samps, int_file=args.int_file, pep_colname=args.pep_colname,
               data_dir=args.data_dir, overwrite=args.overwrite, outfile=args.outfile, func_file=args.func_file,
               func_colname=args.func_colname, ontology=args.ontology, slim_down=args.slim_down, tax_file=args.tax_file,
               tax_colname=args.tax_colname, nopep=args.nopep, nopep_file=args.nopep_file,
               ft_func_data_dir=args.ft_func_data_dir, ft_tax_data_dir=args.ft_tax_data_dir,
               ft_tar_rank=args.ft_tar_rank)
    elif args.command == "filter":
        run_filter(file=args.expand_file, sinfo=args.samps,
                   ontology=args.ontology, mode=args.mode,
                   qthreshold=args.threshold,
                   min_child_non_leaf=args.min_children_non_leaf, min_child_nsamp=args.min_child_nsamp,
                   min_peptides=args.min_peptides, min_pep_nsamp=args.min_pep_nsamp, outfile=args.outfile)
    elif args.command == "test":
        stat(infile=args.file, samps=args.samps, paired=args.paired, parametric=args.parametric,
             ontology=args.ontology, mode=args.mode, outfile=args.outfile)
    elif args.command == "viz":
        print('viz')
    sys.exit(0)


def parse_args_cli():
    parser = argparse.ArgumentParser()

    # split this into three submodules
    subparsers = parser.add_subparsers(title="commands", dest="command")
    parser_expand = subparsers.add_parser('expand')
    parser_filter = subparsers.add_parser('filter')
    parser_stat = subparsers.add_parser('stat')
    parser_viz = subparsers.add_parser('viz')

    # samps file is required in all four parsers
    for par in (parser_expand, parser_filter, parser_stat, parser_viz):
        common_tmp = par.add_argument_group('Arguments common to all modules.')
        # todo: make example of tabular samps file
        common_tmp.add_argument('--samps', '-s', required=True,
                                help='Give the column names in the intensity file that ' +
                                     'correspond to a given sample group. ' +
                                     'This can either be JSON formatted or be a path to a tabular file. ' +
                                     'JSON example of two experimental groups and two samples in each group: ' +
                                     '{"A": ["A1", "A2"], "B": ["B1", "B2"]}')
        common_tmp = par.add_argument_group('Arguments for all analyses')
        common_tmp.add_argument('--mode', '-m', choices=['fn', 'tax', 'taxfn'], required=True,
                            help='Analysis mode. If taxfn is chosen, both function and taxonomy files must be provided')
        common_tmp.add_argument('--ontology', choices=['go', 'cog', 'ec'], required=False,
                          help='Which functional terms to use. Ignored (and not required) if mode is not fn or taxfn.')

    # ---- METAQUANTOME EXPAND ---- #
    common = parser_expand.add_argument_group('Arguments for all 3 modes')
    common.add_argument('--nopep', action="store_true",
                        help="If provided, need to provide a --nopep_file.")
    common.add_argument('--nopep_file',
                        help="File with functional annotations and intensities. ")
    common.add_argument('--int_file', '-i',
                        help='Path to the file with intensity data. Must be tabular, have a peptide sequence column, '+
                             'and be raw, untransformed intensity values. Missing values can be 0, NA, or NaN' +
                             '- transformed to NA for analysis')
    common.add_argument('--pep_colname',
                        help='The column name within the intensity, function, and/or taxonomy file that corresponds ' +
                             'to the peptide sequences. ')
    common.add_argument('--outfile', required=True,
                        help='Output file')

    f_or_t = parser_expand.add_argument_group('Function for function or taxonomy alone')
    f_or_t.add_argument('--data_dir',
                        help='Path to database directory. Pre-downloaded databases can be stored in a separate' +
                             ' directory and timestamped. '+
                             'Note that names of files within the directory cannot be changed. ' +
                             'The default is the mode-appropriate subdirectory of <metaquant_package_root>/data.')


    # function-specific
    func = parser_expand.add_argument_group('Function')
    func.add_argument('--func_file', '-f',
                      help='Path to file with function. The file must be tabular, with a peptide sequence column '+
                           'and either a GO-term column, COG column, or EC number column. The name of the functional'
                           ' column should be given in --func_colname. Other columns will be ignored. ')
    func.add_argument('--func_colname',
                      help='Name of the functional column')
    func.add_argument('--slim_down', action='store_true',
                      help='Flag. If provided, terms are mapped from the full OBO to the slim OBO. ' +
                           'Terms not in the full OBO will be skipped.')
    func.add_argument('--overwrite', action='store_true',
                      help='Flag. The most relevant database (GO or EC) is downloaded to data_dir, ' +
                           'overwriting any previously downloaded databases at these locations.')

    # taxonomy-specific
    tax = parser_expand.add_argument_group('Taxonomy')
    tax.add_argument('--tax_file', '-t',
                     help='Path to (tabular) file with taxonomy assignments. There should be a peptide sequence ' +
                          'column with name pep_colname, and a taxonomy column with name tax_colname')
    tax.add_argument('--tax_colname',
                     help='Name of taxonomy column in tax file. The column ' +
                          'must be either NCBI taxids (strongly preferred) or taxonomy names. ' +
                          'Unipept name output is known to function well, but other formats may not work.')

    # function-taxonomy
    ft = parser_expand.add_argument_group('Function-Taxonomy')
    ft.add_argument('--ft_func_data_dir',
                    help="Path to function data directory for the taxfn mode. " +
                         "The default is <metaquant_pkg_dir>/data/go")
    ft.add_argument('--ft_tax_data_dir',
                    help="Path to function data directory for the taxfn mode.")
    ft.add_argument('--ft_tar_rank', default='genus',
                    help="Desired rank for taxonomy. The default is 'genus'.")

    # ---- METAQUANTOME FILTER ---- #
    parser_filter.add_argument('--expand_file',
                               help="Output from metaquantome expand.")
    parser_filter.add_argument('--min_peptides', default=0, type=int,
                               help='Used for filtering to well-supported annotations. The number of peptides providing ' +
                                    'evidence for a term is the number of peptides directly annotated with that term ' +
                                    'plus the number of peptides annotated with any of its descendants. ' +
                                    'Terms with a number of peptides greater than or equal to min_peptides are retained. ' +
                                    'The default is 0.')
    parser_filter.add_argument('--min_pep_nsamp', default='all',
                               help="Number of samples per group that must meet or exceed min_peptides. " +
                               "Can either be a nonnegative integer or 'all'.")
    parser_filter.add_argument('--min_children_non_leaf', default=0, type=int,
                               help='Used for filtering to informative annotations. ' +
                                    'A term is retained if it has a number of children ' +
                                    'greater than or equal to min_children_non_leaf. ' +
                                    'The default is 0. ')
    parser_filter.add_argument('--min_child_nsamp', default='all',
                               help="Number of samples per group that must meet or exceed min_children_nsamp. " +
                                    "Can either be a nonnegative integer or 'all'.")
    parser_filter.add_argument('--qthreshold', type=int, default=3,
                               help='Minimum number of intensities in each sample group. ' +
                                    'Any functional/taxonomic term with lower number of per-group intensities ' +
                                    'will be filtered out. The default is 3, because this is the minimum ' +
                                    'number for t-tests.')
    parser_filter.add_argument('--outfile', required=True,
                               help="Output file")

    # ---- METAQUANTOME STAT ---- #
    # statistics
    parser_stat.add_argument('--file', '-f', required=True,
                             help='Output file from metaquantome expand.')
    parser_stat.add_argument('--outfile', required=True,
                             help='Output file')
    parser_stat.add_argument('--parametric', type=bool, default=False,
                             help='Choose the type of test. If --parametric True is provided,' +
                                  'then a standard t-test is performed. ' +
                                  'If --parametric False (the default) is provided, ' +
                                  'then a Wilcoxon test is performed.')
    parser_stat.add_argument('--paired', action='store_true',
                             help='Perform paired tests.')

    # ---- METAQUANTOME VIZ ---- #

    # todo

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    cli()
