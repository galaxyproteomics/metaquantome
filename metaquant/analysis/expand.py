from metaquant.util import io
from metaquant.analysis.functional_analysis import functional_analysis
from metaquant.analysis.taxonomy_analysis import taxonomy_analysis
from metaquant.analysis.function_taxonomy_interaction import function_taxonomy_analysis
from metaquant.SampleGroups import SampleGroups


def expand(mode, samps, int_file=None, pep_colname='peptide', data_dir=None, overwrite=False, outfile=None,
           func_file=None, func_colname=None, ontology='go', slim_down=False, tax_file=None, tax_colname=None,
           nopep=False, nopep_file=None, ft_func_data_dir=None, ft_tax_data_dir=None, ft_tar_rank='genus'):

    samp_grps = SampleGroups(samps)

    # read and join files - depending on pep/nopep
    if nopep:
        df = io.read_nopep_table(file=nopep_file, samp_grps=samp_grps,
                                 mode=mode,func_colname=func_colname,
                                 tax_colname=tax_colname)
    else:
        df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                    int_file=int_file,
                                    func_file=func_file, func_colname=func_colname,
                                    tax_colname=tax_colname, tax_file=tax_file)
    # run analysis based on modes
    if mode == 'fn':
        results = functional_analysis(df=df, func_colname=func_colname, samp_grps=samp_grps, ontology=ontology,
                                      slim_down=slim_down, data_dir=data_dir, overwrite=overwrite)
    elif mode == 'tax':
        results = taxonomy_analysis(df=df, samp_grps=samp_grps, data_dir=data_dir, tax_colname=tax_colname)
    elif mode == 'taxfn':
        results = function_taxonomy_analysis(df=df, func_colname=func_colname, pep_colname=pep_colname,
                                             ontology=ontology, overwrite=overwrite, slim_down=slim_down,
                                             tax_colname=tax_colname, samp_grps=samp_grps, ft_tar_rank=ft_tar_rank,
                                             ft_func_data_dir=ft_func_data_dir, ft_tax_data_dir=ft_tax_data_dir)
    else:
        raise ValueError("Invalid mode. Expected one of: %s" % ['fun', 'tax', 'taxfn'])
    # filter min observed here

    # set up written output
    if outfile:
        cols = io.define_outfile_cols_expand(samp_grps, ontology, mode)
        results.to_csv(outfile, columns=cols, sep="\t", header=True, index=False, na_rep="NA")
    # whether writing out or not, return result data frame
    return results
