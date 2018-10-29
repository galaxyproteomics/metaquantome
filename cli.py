import sys
import argparse
from metaquant.analysis.expand import expand


def cli():
    args = parse_args_cli()
    # TODO: pass args to metaquant_runner as a list
    # todo: also maybe split metaquant_runner into three functions to match the cli

    if args.command == "expand":
        print('expand')
        expand(args)
        # run expand function
    elif args.command == "test":
        print('test')
        # run test function
    elif args.command == "viz":
        print('viz')
        # run expand function

    # runner.metaquant_runner(args.mode, sinfo=args.samps, int_file=args.int_file, pep_colname=args.pep_colname,
    #                         func_file=args.func_file, func_colname=args.func_colname, tax_file=args.tax_file,
    #                         ontology=args.ontology, tax_colname=args.tax_colname, outfile=args.outfile,
    #                         slim_down=args.slim_down, test=args.test,
    #                         paired=args.paired, threshold=args.threshold, data_dir=args.data_dir, overwrite=args.overwrite,
    #                         min_peptides=args.min_peptides, min_children_non_leaf=args.min_children_non_leaf)
    sys.exit(0)


def parse_args_cli():
    parser = argparse.ArgumentParser()

    # split this into three submodules
    subparsers = parser.add_subparsers(title="commands", dest="command")
    parser_expand = subparsers.add_parser('expand')
    parser_test = subparsers.add_parser('test')
    parser_viz = subparsers.add_parser('viz')

    # we need these two arguments in all three parsers
    for par in (parser_expand, parser_test, parser_viz):
        common_tmp = par.add_argument_group('Required for all analyses')
        common_tmp.add_argument('--mode', '-m', choices=['fn', 'tax', 'taxfn'], required=True,
                            help='Analysis mode. If taxfun is chosen, both function and taxonomy files must be provided')
        # todo: make example of tabular samps file
        common_tmp.add_argument('--samps', '-s', required=True,
                            help='Give the column names in the intensity file that ' +
                                 'correspond to a given sample group. ' +
                                 'This can either be JSON formatted or be a path to a tabular file. ' +
                                 'JSON example of two experimental groups and two samples in each group: ' +
                                 '{"A": ["A1", "A2"], "B": ["B1", "B2"]}')

    # METAQUANTOME EXPAND #
    common = parser_expand.add_argument_group('Arguments for all 3 modes')
    common.add_argument('--pep', type=bool, default=True,
                        help="Peptides or no peptides.")
    common.add_argument('--int_file', '-i', required=True,
                        help='Path to the file with intensity data. Must be tabular, have a peptide sequence column, '+
                             'and be raw, untransformed intensity values. Missing values can be 0, NA, or NaN' +
                             '- transformed to NA for analysis')
    common.add_argument('--pep_colname', required=True,
                        help='The column name within the intensity, function, and/or taxonomy file that corresponds ' +
                             'to the peptide sequences. ')
    common.add_argument('--outfile', required=True,
                        help='Output file')
    common.add_argument('--data_dir',
                        help='Path to database directory. Pre-downloaded databases can be stored in a separate' +
                             ' directory and timestamped. '+
                             'Note that names of files within the directory cannot be changed. ' +
                             'The default is <metaquant_package_root>/data.')
    common.add_argument('--min_peptides', default=0, type=int,
                        help='Used for filtering to well-supported annotations. The number of peptides providing ' +
                             'evidence for a term is the number of peptides directly annotated with that term ' +
                             'plus the number of peptides annotated with any of its descendants. ' +
                             'Terms with a number of peptides greater than or equal to min_peptides are retained.')
    common.add_argument('--min_children_non_leaf', default=0, type=int,
                        help='Used for filtering to informative annotations. ' +
                             'A term is retained if it has a number of children ' +
                             'greater than or equal to min_children_non_leaf. ')
    common.add_argument('--threshold', type=int, default=3,
                        help='Minimum number of intensities in each sample group. ' +
                             'Any functional/taxonomic term with lower number of per-group intensities ' +
                             'will be filtered out. The default is 3, because this is the minimum ' +
                             'number for t-tests.')

    # function-specific
    func = parser_expand.add_argument_group('Function')
    func.add_argument('--func_file', '-f',
                      help='Path to file with function. The file must be tabular, with a peptide sequence column '+
                           'and either a GO-term column, COG column, or EC number column. The name of the functional'
                           ' column should be given in --func_colname. Other columns will be ignored. ')
    func.add_argument('--func_colname',
                      help='Name of the functional column')
    func.add_argument('--ontology', choices=['go', 'cog', 'ec'],
                      help='Which functional terms to use.')
    func.add_argument('--slim_down', action='store_true',
                      help='Flag. If provided, terms are mapped from the full OBO to the slim OBO. ' +
                           'Terms not in the full OBO will be skipped.')
    func.add_argument('--overwrite', action='store_true',
                      help='Flag. If provided, the most relevant databases (GO and/or EC) are downloaded to data_dir, '+
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

    # METAQUANTOME TEST #

    # statistics
    parser_test.add_argument('--file', '-f', help='Output file from metaquantome expand.')
    parser_test.add_argument('--parametric', type=bool, default=True,
                      help='Choose the type of test. If --parametric True is provided,' +
                           'then a standard t-test is performed. If --parametric False is provided, ' +
                           'then a Wilcoxon test is performed.')
    parser_test.add_argument('--paired', action='store_true',
                      help='Perform paired tests.')

    # METAQUANTOME VIZ #

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    cli()
