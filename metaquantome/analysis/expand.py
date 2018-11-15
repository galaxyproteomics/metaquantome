import numpy as np

import metaquantome.util.expand_io as expand_io
from metaquantome.SampleAnnotations import SampleAnnotations
from metaquantome.analysis.functional_analysis import functional_analysis
from metaquantome.analysis.taxonomy_analysis import taxonomy_analysis
from metaquantome.analysis.function_taxonomy_interaction import function_taxonomy_analysis
from metaquantome.SampleGroups import SampleGroups


def expand(mode, samps, int_file=None, pep_colname='peptide', data_dir=None, overwrite=False, outfile=None,
           func_file=None, func_colname=None, ontology='go', slim_down=False, tax_file=None, tax_colname=None,
           nopep=False, nopep_file=None, ft_func_data_dir=None, ft_tax_data_dir=None, ft_tar_rank='genus'):
    # todo: doc
    samp_grps = SampleGroups(samps)

    # read and join files - depending on pep/nopep
    if nopep:
        df = expand_io.read_nopep_table(file=nopep_file, samp_grps=samp_grps,
                                        mode=mode, func_colname=func_colname,
                                        tax_colname=tax_colname)
    else:
        df = expand_io.read_and_join_files(mode, pep_colname, samp_grps,
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
        cols = expand_io.define_outfile_cols_expand(samp_grps, ontology, mode)
        expand_io.write_out_general(results, outfile=outfile, cols=cols)
    # whether writing out or not, return result data frame
    return results


def common_hierarchical_analysis(db, df, annot_colname, samp_grps, hierarchical=True):

    # import pdb; pdb.set_trace()
    if hierarchical:
        samp_annot = SampleAnnotations(db)
        # make a hierarchy for each sample
        samp_annot.add_samples_from_df(df, annot_colname, samp_grps)
        intensity_all_ranks = samp_annot.to_dataframe()
    else:
        intensity_all_ranks = df

    # calculate means
    int_w_means = calc_means(intensity_all_ranks, samp_grps)
    # todo: calculate sds

    # clean and log transform
    # replace nan with zero, so that np.log2 returns nan
    int_w_means[int_w_means == 0] = np.nan

    # take log of intensities for return
    int_w_means[samp_grps.all_intcols] = np.log2(int_w_means[samp_grps.all_intcols])

    int_w_means['id'] = int_w_means.index

    return int_w_means


def calc_means(df, samp_grps):

    for i in range(samp_grps.ngrps):
        grp_name = samp_grps.grp_names[i]
        samples_in_grp = samp_grps.sample_names[grp_name]
        if len(samples_in_grp) > 1:
            sample = df[samples_in_grp]
            means = np.mean(sample, axis = 1)
            means[means == 0] = np.nan
            logs = np.log2(means)
            df[samp_grps.mean_names[i]] = logs

        else:
            # just log transform the single sample
            df[samp_grps.mean_names[i]] = np.log2(df[samples_in_grp])

    return df