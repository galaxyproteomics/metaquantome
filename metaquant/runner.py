from metaquant.SampleGroups import SampleGroups
from metaquant.analysis.function_taxonomy_interaction import function_taxonomy_analysis
from metaquant.analysis.functional_analysis import functional_analysis
from metaquant.analysis.taxonomy_analysis import taxonomy_analysis
from metaquant.util import io
import os
import logging
import sys


def runner(mode, sinfo, int_file, pep_colname='peptide', func_colname=None, func_file=None, tax_file=None,
           ontology='go', tax_colname=None, outfile=None, slim_down=False, test=False, paired=False, parametric=False,
           threshold=0, data_dir=None, overwrite=False, min_peptides=0, min_children_non_leaf=0):

    # initialize logger
    logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stderr)

    # define object with sample groups, intensity columns, etc.
    samp_grps = SampleGroups(sinfo)
    if not os.path.exists(int_file):
        raise FileNotFoundError('int_file does not exist. Please check filename')

    # make sure base data directory exists, if provided
    if data_dir:
        if not os.path.isdir(data_dir):
            print('data_dir is not a directory or was not found. Making data_dir')
            os.mkdir(data_dir)

    df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                int_file=int_file, func_file=func_file, func_colname=func_colname,
                                tax_colname=tax_colname, tax_file=tax_file)

    # run analysis based on modes
    if mode == 'fn':
        results = functional_analysis(df=df, func_colname=func_colname, samp_grps=samp_grps, test=test,
                                      threshold=threshold, ontology=ontology, slim_down=slim_down, paired=paired,
                                      parametric=parametric, data_dir=data_dir, overwrite=overwrite,
                                      min_peptides=min_peptides, min_children_non_leaf=min_children_non_leaf)
    elif mode == 'tax':
        results = taxonomy_analysis(df=df, samp_grps=samp_grps, test=test, threshold=threshold, paired=paired,
                                    parametric=parametric, data_dir=data_dir, tax_colname=tax_colname,
                                    min_peptides=min_peptides, min_children_non_leaf=min_children_non_leaf)
    elif mode == 'taxfn':
        if ontology != 'cog':
            raise ValueError("Only cog is supported for ft interaction. " +
                             "Make sure you have a cog column and supply the column name to func_colname")
        results = function_taxonomy_analysis(df=df, cog_name=func_colname, lca_colname=tax_colname, samp_grps=samp_grps,
                                             test=test, threshold=threshold, paired=paired, parametric=parametric,
                                             data_dir=data_dir)
    else:
        raise ValueError("Invalid mode. Expected one of: %s" % ['fun', 'tax', 'taxfn'])
    # set up written output
    if outfile:
        cols = define_outfile_cols(samp_grps, ontology, mode, test)
        results.to_csv(outfile, columns=cols, sep="\t", header=True, index=False, na_rep="NA")
    # whether writing out or not, return result data frame
    return results


def define_outfile_cols(samp_grps, ontology, mode, test):
    int_cols = []
    if test:
        fc_name = ['log2fc_' + samp_grps.grp_names[0] + '_over_' + samp_grps.grp_names[1]]
        int_cols += fc_name + ['p', 'corrected_p']
    int_cols += samp_grps.mean_names + samp_grps.all_intcols
    if mode == 'fn':
        if ontology == 'go':
            cols = ['id', 'name', 'namespace'] + int_cols
        elif ontology == 'cog':
            cols = ['id', 'description'] + int_cols
        elif ontology == 'ec':
            cols = ['id', 'description'] + int_cols
        else:
            raise ValueError("Invalid ontology. Expected one of: %s" % ['go', 'cog', 'ec'])
    elif mode == 'tax':
        cols = ['id', 'taxon_name', 'rank'] + int_cols
    elif mode == 'taxfn':
        cols = [ontology, 'cog_descript', 'taxon_name', 'rank'] + int_cols
    else:
        raise ValueError("Invalid mode. Expected one of: %s" % ['fun', 'tax', 'taxfn'])
    return cols

