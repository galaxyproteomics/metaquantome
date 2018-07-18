from metaquant.function_taxonomy_interaction import function_taxonomy_analysis
from metaquant.functional_analysis import functional_analysis
from metaquant.taxonomy_analysis import taxonomy_analysis
from metaquant import io
from metaquant.definitions import DATA_DIR
import os
import logging
import sys


def metaquant(mode, sample_names, int_file, pep_colname='peptide', func_colname=None, func_file=None, tax_file=None,
              ontology='go', tax_colname=None, outfile=None, slim_down=False, test=False, paired=False, threshold=0,
              data_dir=None, overwrite=False):

    # initialize logger
    logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stderr)

    # define base data directory
    if not data_dir:
        data_dir = DATA_DIR
    else:
        if not os.path.isdir(data_dir):
            raise NotADirectoryError('data_dir is not a directory or was not found.')

    # define object with sample groups, intensity columns, etc.
    samp_grps = io.SampleGroups(sample_names)

    if not os.path.exists(int_file):
        raise FileNotFoundError('int_file does not exist. Please check filename')

    modes = ['fn', 'tax', 'taxfn']
    if mode == 'fn':
        function_check(func_file, func_colname)
        df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                    int_file=int_file, func_file=func_file, func_colname=func_colname)
        results = functional_analysis(df=df, func_colname=func_colname, samp_grps=samp_grps, test=test, threshold=threshold,
                                      ontology=ontology, slim_down=slim_down, paired=paired, data_dir=data_dir,
                                      overwrite=overwrite)

    elif mode == 'tax':
        tax_check(tax_file, tax_colname)
        df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                    int_file=int_file, tax_file=tax_file, tax_colname=tax_colname)
        results = taxonomy_analysis(df=df, samp_grps=samp_grps, test=test, threshold=threshold, paired=paired,
                                    data_dir=data_dir, tax_colname=tax_colname)

    elif mode == 'taxfn':
        function_check(func_file, func_colname)
        tax_check(tax_file, tax_colname)
        if ontology != 'cog':
            raise ValueError("Only cog is supported for ft interaction. " +
                             "Make sure you have a cog column and supply the column name to func_colname")

        df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                    int_file=int_file,
                                    tax_file=tax_file, tax_colname=tax_colname,
                                    func_file=func_file, func_colname=func_colname)
        results = function_taxonomy_analysis(df=df, cog_name=func_colname, lca_colname=tax_colname, samp_grps=samp_grps,
                                             test=test, threshold=threshold, paired=paired, data_dir=data_dir,
                                             overwrite=overwrite)

    else:
        raise ValueError("Invalid mode. Expected one of: %s" % modes)

    if outfile:
        cols = []
        int_cols = []

        if test:
            fc_name = ['log2fc_' + samp_grps.grp_names[0] + '_over_' + samp_grps.grp_names[1]]
            int_cols += fc_name + ['p', 'corrected_p']

        int_cols += samp_grps.mean_names + samp_grps.all_intcols

        if mode == 'fn':
            if ontology == 'go':
                cols = ['go_id', 'name', 'namespace'] + int_cols
            if ontology == 'cog':
                cols = ['cog', 'description'] + int_cols

        if mode == 'tax':
            cols = ['taxon_name', 'rank'] + int_cols

        if mode == 'taxfn':
            cols = [ontology, 'cog_descript', 'taxon_name', 'rank'] + int_cols

        results.to_csv(outfile, columns=cols, sep="\t", header=True, index=False, na_rep="NA")

    return results


def function_check(func_file, func_colname):
    if not func_file:
        raise IOError('Function tabular file not provided (--func_file)')
    if not os.path.exists(func_file):
        raise FileNotFoundError('func_file does not exist. Please check filename')
    if not func_colname:
        raise ValueError('func_colname=None. Please provide a function column name (--func_colname)')


def tax_check(tax_file, tax_colname):
    if not tax_file:
        raise IOError('Taxonomy tabular file not provided (--tax_file)')
    if not os.path.exists(tax_file):
        raise FileNotFoundError('func_file does not exist. Please check filename')
    if not tax_colname:
        raise ValueError('tax_colname=None. Please provide a taxonomy column name (--tax_colname)')
