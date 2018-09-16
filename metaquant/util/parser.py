import argparse


def parse_args_cli():
    parser = argparse.ArgumentParser()

    # required for all analyses
    common = parser.add_argument_group('Required for all analyses')
    common.add_argument('--mode', '-m', choices=['fn', 'tax', 'taxfn'], required=True,
                        help='Analysis mode. If taxfun is chosen, both function and taxonomy files must be provided')
    common.add_argument('--samps', '-s', required=True,
                        help='Give the column names in the intensity file that ' +
                             'correspond to a given sample group. ' +
                             'This can either be JSON formatted or be a path to a tabular file. ' +
                             'JSON example of two experimental groups and two samples in each group: ' +
                             '{"A": ["A1", "A2"], "B": ["B1", "B2"]}')
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

    # function-specific
    func = parser.add_argument_group('Function')
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
    tax = parser.add_argument_group('Taxonomy')
    tax.add_argument('--tax_file', '-t',
                     help='Path to (tabular) file with taxonomy assignments. There should be a peptide sequence ' +
                          'column with name pep_colname, and a taxonomy column with name tax_colname')
    tax.add_argument('--tax_colname',
                     help='Name of taxonomy column in tax file. The column ' +
                          'must be either NCBI taxids (strongly preferred) or taxonomy names. ' +
                          'Unipept name output is known to function well, but other formats may not work.')

    # statistics
    stat = parser.add_argument_group('Statistics')
    stat.add_argument('--test', action='store_true',
                      help='Perform t-tests on the summed intensities.' +
                           'FDR-corrected q-values are returned.')
    stat.add_argument('--paired', action='store_true',
                      help='If --test and --paired are provided, perform paired t-tests.')
    stat.add_argument('--threshold', type=int,
                      help='Minimum number of intensities in each sample group. ' +
                           'Anything with lower number of intensities will be filtered out.' +
                           'Only applies when testing, not to descriptive statistics.')

    args = parser.parse_args()
    return args
