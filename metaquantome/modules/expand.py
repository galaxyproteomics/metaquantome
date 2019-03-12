import numpy as np

import metaquantome.util.expand_io as expand_io
from metaquantome.classes.SampleAnnotations import SampleAnnotations
from metaquantome.modules.functional_analysis import functional_analysis
from metaquantome.modules.taxonomy_analysis import taxonomy_analysis
from metaquantome.modules.function_taxonomy_interaction import function_taxonomy_analysis
from metaquantome.classes.SampleGroups import SampleGroups


def expand(mode, sinfo, int_file, pep_colname_int, pep_colname_func, pep_colname_tax, data_dir=None, outfile=None,
           func_file=None, func_colname=None, ontology='go', slim_down=False, tax_file=None, tax_colname=None,
           nopep=False, nopep_file=None, ft_tar_rank='genus'):
    """
    Expand the directly annotated hierarchy to one with all ancestors.

    :param mode: either 't', 'f', or 'ft'
    :param sinfo: Either a json string with experimental information or a path to a tabular file
    :param int_file: Path to the tabular intensity file
    :param pep_colname_int: peptide column name in int_file
    :param pep_colname_func: peptide column name in func_file
    :param pep_colname_tax: peptide column name in tax_file
    :param outfile: path to write results to
    :param func_file: path to functional annotations
    :param func_colname: column name for the functional annotations within func_file
    :param ontology: for function mode only. either 'go', 'ec', or 'cog'
    :param slim_down: if True, maps full GO terms to slim
    :param tax_file: path to taxonomy file
    :param tax_colname: column name with taxonomy annotations
    :param nopep: if True, do nopep modules
    :param nopep_file: path to file without peptides
    :param data_dir: path to parent directory of database files
    :param ft_tar_rank: in ft mode, all taxonomy are mapped to this rank if possible.
    :return: returns a dataframe of functional or taxonomic terms with associated intensities.
    Missing values are represented as 0.
    """
    # define the sample groups object
    samp_grps = SampleGroups(sinfo)

    # read and join files - depending on pep/nopep
    if nopep:
        df = expand_io.read_nopep_table(file=nopep_file, samp_grps=samp_grps,
                                        mode=mode, func_colname=func_colname,
                                        tax_colname=tax_colname)
    else:
        df = expand_io.read_and_join_files(mode, pep_colname_int=pep_colname_int, pep_colname_func=pep_colname_func,
                                           pep_colname_tax=pep_colname_tax, samp_grps=samp_grps, int_file=int_file,
                                           tax_file=tax_file, func_file=func_file, func_colname=func_colname,
                                           tax_colname=tax_colname)
    # run modules based on modes
    if mode == 'f':
        results = functional_analysis(df=df, func_colname=func_colname, samp_grps=samp_grps, ontology=ontology,
                                      slim_down=slim_down, data_dir=data_dir)
    elif mode == 't':
        results = taxonomy_analysis(df=df, samp_grps=samp_grps, data_dir=data_dir, tax_colname=tax_colname)
    elif mode == 'ft':
        results = function_taxonomy_analysis(df=df, func_colname=func_colname, pep_colname=pep_colname_int,
                                             ontology=ontology, slim_down=slim_down, tax_colname=tax_colname,
                                             samp_grps=samp_grps, ft_tar_rank=ft_tar_rank, data_dir=data_dir)
    else:
        raise ValueError('Invalid mode. Expected one of "f", "t", or "ft"')

    # set up written output
    if outfile:
        cols = expand_io.define_outfile_cols_expand(samp_grps, ontology, mode)
        expand_io.write_out_general(results, outfile=outfile, cols=cols)
    # whether writing out or not, return result data frame (mostly for testing)
    return results


def common_hierarchical_analysis(db, df, annot_colname, samp_grps, hierarchical=True):
    """
    Create hierarchies from original dataframe

    :param db: reference database
    :param df: original dataframe.
    :param annot_colname: column name containing either the functional or taxonomic annotations
    :param samp_grps: SampleGroups object
    :param hierarchical: False should be used in the case of COG, in which case it just adds up the intensities
    for each term
    :return: DataFrame with summarised intensities and other quantitative measures
    """
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

    # define an "id" column using the index
    int_w_means['id'] = int_w_means.index
    return int_w_means


def calc_means(df, samp_grps):
    """
    Calculate the groupwise average for each term
    :param df: expanded dataframe, non-transformed intensities
    :param samp_grps: SampleGroups() object
    :return: dataframe with a mean column, which is the log of the mean of the group-specific intensities
    """
    for i in range(samp_grps.ngrps):
        grp_name = samp_grps.grp_names[i]
        samples_in_grp = samp_grps.sample_names[grp_name]
        if len(samples_in_grp) > 1:
            sample = df[samples_in_grp]
            means = np.mean(sample, axis=1)
            means[means == 0] = np.nan
            logs = np.log2(means)
            df[samp_grps.mean_names[i]] = logs

        else:
            # just log transform the single sample
            df[samp_grps.mean_names[i]] = np.log2(df[samples_in_grp])

    return df
