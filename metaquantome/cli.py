import sys
import argparse
import logging
import os
import pkg_resources

# add metaquantome parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from metaquantome.modules.expand import expand
from metaquantome.modules.filter import run_filter
from metaquantome.modules.stat import stat
from metaquantome.modules.run_viz import run_viz
from metaquantome.modules.db_download_handler import db_download_handler


def cli():
    """
    Command line interface; main entry point to metaQuantome

    :return: exit code
    """
    # initialize logger
    logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stderr)
    args = parse_args_cli()
    if args.command == "db":
        db_download_handler(args.dbs, args.dir, args.update)
    elif args.command == "expand":
        expand(mode=args.mode, sinfo=args.samps, int_file=args.int_file, pep_colname_int=args.pep_colname_int,
               pep_colname_func=args.pep_colname_func, pep_colname_tax=args.pep_colname_tax, data_dir=args.data_dir,
               outfile=args.outfile, func_file=args.func_file, func_colname=args.func_colname, ontology=args.ontology,
               slim_down=args.slim_down, tax_file=args.tax_file, tax_colname=args.tax_colname, nopep=args.nopep,
               nopep_file=args.nopep_file, ft_tar_rank=args.ft_tar_rank)
    elif args.command == "filter":
        run_filter(expanded_file=args.expand_file, sinfo=args.samps, ontology=args.ontology, mode=args.mode,
                   qthreshold=args.qthreshold, min_child_non_leaf=args.min_children_non_leaf,
                   min_child_nsamp=args.min_child_nsamp, min_peptides=args.min_peptides,
                   min_pep_nsamp=args.min_pep_nsamp, outfile=args.outfile)
    elif args.command == "stat":
        stat(infile=args.file, sinfo=args.samps, paired=args.paired, parametric=args.parametric, ontology=args.ontology,
             mode=args.mode, outfile=args.outfile, control_group=args.control_group)
    elif args.command == "viz":
        run_viz(plottype=args.plottype,
                img=args.img,
                infile=args.infile,
                mode=args.mode,
                meancol=args.meancol,
                nterms=args.nterms,
                strip=args.strip,
                target_rank=args.target_rank,
                barcol=args.barcol,
                textannot=args.textannot,
                calculate_sep=args.calculate_sep,
                fc_name=args.fc_name,
                fc_corr_p=args.fc_corr_p,
                flip_fc=args.flip_fc,
                gosplit=args.gosplit,
                sinfo=args.samps,
                filter_to_sig=args.filter_to_sig,
                alpha=args.alpha,
                whichway=args.whichway,
                name=args.name,
                id=args.id,
                target_onto=args.target_onto,
                width=args.width,
                height=args.height,
                tabfile=args.tabfile,
                feature_cluster_size=args.feature_cluster_size,
                sample_cluster_size=args.sample_cluster_size)
    else:
        ValueError('incorrect mode. please provide one of "db", "expand", "filter", "stat", or "viz".')
    sys.exit(0)


def check_col_range(arg):
    try:
        value = int(arg)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err))
    if value not in [1, 2, 3, 4, 5, 6]:
        message = "Expected value in [1, 2, 3, 4, 5, 6], got value = {}".format(value)
        raise argparse.ArgumentTypeError(message)
    return value


def parse_args_cli():
    """
    parse the command line arguments
    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="metaQuantome is a tool that performs quantitative analysis"
            " on the function and taxonomy of microbomes and their interactions. "
            "For more background information, please read the associated manuscript: https://doi.org/10.1074/mcp.ra118.001240. "
            "For a more hands-on tutorial, please visit the following page: https://galaxyproteomics.github.io/metaquantome_mcp_analysis/cli_tutorial/cli_tutorial.html.\n\n"
            "The metaQuantome workflow is as follows: db → expand → filter → stat → viz.\n\n" 
            "Run `metaquantome {db,expand,filter,stat,viz} -h` for more information on the individual modules. "
            "Any issues can be brought to attention here: https://github.com/galaxyproteomics/metaquantome/issues."
    )
    # parser.add_argument('--help', action=_HelpAction, help="Help")
    parser.add_argument('-v', '--version', action='version',
                         version=pkg_resources.require("metaquantome")[0].version)

    # split this into three submodules
    subparsers = parser.add_subparsers(title="commands", dest="command")
    parser_db = subparsers.add_parser('db',
        formatter_class=argparse.RawTextHelpFormatter,
        description="metaQuantome uses freely available bioinformatic databases to expand your set of direct annotations. For most cases, all 3 databases can be downloaded (the default).\n"
                    "The databases are:\n"
                    "1. NCBI taxonomy database. This contains a list of all currently identified taxa and the relationships between them.\n"
                    "2. Gene Ontology (GO) term database. metaQuantome uses the OBO format of the database. Specifically, two files are used: the go-basic.obo file, which is a simplified version of the GO database that is guaranteed to be acyclic, and the metagenomics slim GO, which is a subset of the full GO that is useful for microbiome research. More details are available at http://geneontology.org/docs/download-ontology/\n"
                    "3. ENZYME database with Enzyme Classification (EC) numbers. This database classifies enzymes and organizes the relationships between them.\n"
                    "This module downloads the most recent releases of the specified databases and stores them in a single file, which can then be accessed by the rest of the metaQuantome modules. For reference, the taxonomy database is the largest (~500 Mb), while the GO and EC databases are smaller: ~34 Mb and ~10Mb, respectively."
                    )
    parser_expand = subparsers.add_parser('expand', 
        formatter_class=argparse.RawTextHelpFormatter,
        description='The expand module is the first analysis step in the metaQuantome analysis workflow, and can be run to analyze function, taxonomy, or function and taxonomy together.\n'
                    'To prepare to run this module, you must create your samples file (tabular file with `group` and `colnames` columns with each row as a sample group) and download the necessary databases with `metaquantome db`.\n'
                    'Some example analysis workflows are:\n'
                    '1. Get the functional, taxonomic, or functional-taxonomic distribution: run expand, filter, and viz.\n'
                    '2. Cluster analysis: run expand, filter, and viz. The viz module has heatmaps and PCA plots for cluster analysis.\n'
                    '3. Differential expression: run expand, filter, stat, and viz.\n\n'
                    'The following information is required for all 3 analysis modes (function, taxonomy, and function-taxonomy).\n'
                    '- experimental design information.\n'
                    '- a tab-separated peptide intensity file.\n'
                    '- the name of the peptide column in the intensity file.\n\n'
                    'Function mode\n'
                    'In function mode, the following information is required:\n'
                    '- the ontology being used: Gene Ontology (GO), Clusters of Orthologous Groups (COG), or Enzyme Commission (EC) numbers.\n'
                    '- a tab-separated functional annotation file, with a peptide column and a functional annotation column. An entry in the functional annotation column may contain multiple functional annotations separated by commas.\n'
                    '- the name of the peptide column in the functional annotation file.\n'
                    '- the name of the functional annotation column in the functional annotation file.\n\n'
                    'Taxonomy mode\n'
                    'In taxonomy mode, the following information is required:\n'
                    '- a tab-separated taxonomy annotation file, with a peptide column and a taxonomy annotation column. The taxonomic annotations should be the lowest common ancestor (LCA) for each peptide, preferably given as NCBI taxonomy IDs.\n'
                    '- the name of the peptide column in the taxonomic annotation file.\n'
                    '- the name of the taxonomy annotation column in the taxonomy annotation file.\n\n'
                    'Function-Taxonomy mode\n'
                    'In the combined mode, all of the above must be provided. In addition, the "target rank" must be provided, which is the desired taxonomic rank at which to summarize the function/taxonomy results.')
    parser_filter = subparsers.add_parser('filter',
        formatter_class=argparse.RawTextHelpFormatter, 
        description="The filter module is the second step in the metaQuantome analysis workflow. The purpose of the filter module is to filter expanded terms to those that are representative and well-supported by the data. "
                "Please see the manuscript (https://doi.org/10.1074/mcp.ra118.001240) for further details about filtering.")
    parser_stat = subparsers.add_parser('stat', 
        formatter_class=argparse.RawTextHelpFormatter, 
        description="The stat module is the third step in the metaQuantome analysis workflow."
                    " The purpose of the stat module is to perform differential expression analysis between 2 experimental conditions."
                    " metaQuantome offers paired and unpaired tests, as well as parametric and non-parametric options.")
    parser_viz = subparsers.add_parser('viz', 
        formatter_class=argparse.RawTextHelpFormatter, 
        description="The viz module is the final step in the metaQuantome analysis workflow."
            " The available visualizations are:\n"
                "-bar plot\n-volcano plot\n-heatmap\n-PCA plot\n"
                "Please consult the manuscript (https://doi.org/10.1074/mcp.ra118.001240.) for full details on each of these plots.")

    # db download
    parser_db.add_argument('dbs', type=str, nargs='+', choices=['ncbi', 'go', 'ec'],
                           help='database to download. note that COG mode does not require a download due to ' +
                           'its simplicity.')
    parser_db.add_argument('--dir', '-d', type=str,
                           help='data directory for files.')
    parser_db.add_argument('--update', '-u', action="store_true",
                           help='overwrite existing databases if present.')

    # samps file is required in all four non-db parsers
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
        common_tmp.add_argument('--mode', '-m', choices=['f', 't', 'ft'], required=True,
                                help='Analysis mode. If ft is chosen, both function and taxonomy files must be provided')
        common_tmp.add_argument('--ontology', choices=['go', 'cog', 'ec'], required=False,
                                help='Which functional terms to use. Ignored (and not required) if mode is not f or ft.')
        common_tmp.add_argument('--data_dir',
                                help="Path to data directory. The default is <metaquantome_pkg_dir>/data")

    # ---- METAQUANTOME EXPAND ---- #
    common = parser_expand.add_argument_group('Arguments for all 3 modes')
    common.add_argument('--nopep', action="store_true",
                        help="If provided, need to provide a --nopep_file.")
    common.add_argument('--nopep_file',
                        help="File with functional or taxonomic annotations and intensities. ")
    common.add_argument('--int_file', '-i',
                        help='Path to the file with intensity data. Must be tabular, have a peptide sequence column, '+
                             'and be raw, untransformed intensity values. Missing values can be 0, NA, or NaN' +
                             '- transformed to NA for modules')
    common.add_argument('--pep_colname_int',
                        help='The column name within the intensity file that corresponds ' +
                             'to the peptide sequences. ')
    common.add_argument('--pep_colname_func',
                        help='The column name within the function file that corresponds ' +
                             'to the peptide sequences. ')
    common.add_argument('--pep_colname_tax',
                        help='The column name within the taxonomy file that corresponds ' +
                             'to the peptide sequences. ')
    common.add_argument('--outfile', required=True,
                        help='Output file')

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
    parser_stat.add_argument('--control_group', required=True,
                             help='Sample group name of control samples (will be used as denominator for fold change).')
    
    # ---- METAQUANTOME VIZ ---- #
    parser_viz.add_argument('--plottype', '-p', required=True, choices=['bar', 'volcano', 'heatmap', 'pca', 'ft_dist', 'stacked_bar'],
                            help="Select the type of plot to generate.")
    parser_viz.add_argument('--img', required=True,
                            help='Path to the PNG image file (must end in ".png").')
    parser_viz.add_argument('--width', default="5",
                            help="Width of the image in inches. Defaults vary by plot type.")
    parser_viz.add_argument('--height', default="5",
                            help="Height of the image in inches. Defaults vary by plot type.")
    parser_viz.add_argument('--infile', '-i', required=True,
                            help="Input file from stat or filter.")
    parser_viz.add_argument('--strip', default=None,
                            help="Text to remove from column names for plotting.")
    parser_viz.add_argument('--tabfile', default=None,
                            help="Optional. File to write plot table to.")
    parser_viz.add_argument('--fc_corr_p', default=None,
                            help="Name of the corrected p-value column in the stat dataframe. Used while generating volcano plot and while using filter_to_sig in heatmap")
                      
    bar = parser_viz.add_argument_group('Arguments for barplots - including total taxonomy peptide intensity ("bar"), function-taxonomy ' +
                                        'interaction distributions ("ft_dist"), and stacked taxonomy bar plots ("stacked_bar")')
    bar.add_argument('--meancol',
                     help="(Tax bar and FT dist). Mean intensity column name for desired experimental conditio.")
    bar.add_argument('--nterms', default='5',
                     help="(Tax bar, FT dist, and stacked bar). Number of taxa or functional terms to display. The default is 5.")
    bar.add_argument('--barcol', type=check_col_range, default="6",
                     help="(Tax bar and FT dist). Color for the bar fill. The color vector in R is " +
                          'c("dodgerblue", "darkorange", "yellow2", "red2", "darkviolet", "black"), ' +
                          ' so providing a 1 will give the "dodgerblue" color. These same colors are also used in the ' +
                          ' heatmap and PCA plot, so the colors can be tweaked to match. ')
    bar.add_argument('--target_rank',
                     help="(Tax bar, FT dist, and stacked bar). Taxonomic rank to restrict to in the plot. ")
    bar.add_argument('--target_onto', choices=["mf", "bp", "cc"],
                     help="(Function and FT dist bar only) " +
                          "Ontology to restrict to, for function distribution.")
    bar.add_argument("--whichway", choices=["f_dist", "t_dist"],
                     help="(FT dist only) " +
                          "Which distribution - functional distribution for a taxon (f_dist) or " +
                          "taxonomic distribution for a function (t_dist)?")
    bar.add_argument("--name",
                     help="(FT dist only) " +
                          "Provide either a taxonomic or functional term name. Either provide this or an --id.")
    bar.add_argument("--id",
                     help="(FT dist bar only) " +
                          "Taxonomic or functional term id - either a NCBI taxID or a GO term id (GO:XXXXXXX)")

    volc = parser_viz.add_argument_group('Volcano Plot')
    volc.add_argument('--fc_name',
                      help="Name of the fold change column in the stat dataframe.")
    volc.add_argument('--textannot',
                      help="Name of the text annotation column to optionally include in the volcano." +
                           " If missing, no text will be plotted. ")
    volc.add_argument('--gosplit', action="store_true",
                      help="If using GO terms, whether to make one plot for each of BP, CC, and MF.")
    volc.add_argument('--flip_fc', action="store_true",
                      help="Flag. Whether to flip the fold change (i.e., multiply log fold change by -1)")

    heat = parser_viz.add_argument_group('Heatmap')
    heat.add_argument("--filter_to_sig", action="store_true",
                      help="Flag. Only plot significant terms? Necessitates use of results from `test`.")
    heat.add_argument('--alpha', default='0.05',
                      help="If filter_to_sig, the q-value significance level.")
    heat.add_argument('--feature_cluster_size', default='2',
                      help="Number of clusters 'k' to cut the feature dendrogram tree. Default = 2")
    heat.add_argument('--sample_cluster_size', default='2',
                      help="Number of clusters 'k' to cut the sample dendrogram tree. Default = 2")

    pca = parser_viz.add_argument_group('Principal Components Analysis')
    pca.add_argument("--calculate_sep", action="store_true",
                     help="Flag. Calculate separation between groups and include in title?")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    cli()
