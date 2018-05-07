from src.function_taxonomy_interaction import function_taxonomy_interaction_analysis
from src.functional_analysis import functional_analysis
from src.taxonomy_analysis import taxonomy_analysis
from src import common

# TODO
# write arg parser


def metaquant(mode, file, sample1_colnames, sample2_colnames=None,
              cog_colname="cog",
              go_colname="go",
              tax_rank="genus",
              pep_colname="peptide",
              outfile=None, ontology="GO",
              slim_down=False, test=False,
              test_type="modt",
              paired=False, threshold=0,
              obo_path=None, slim_path=None,
              download_obo=False, overwrite_obo=False):

    modes = ['fn', 'tax', 'taxfn', 'fnpred']
    if mode not in modes:
        raise ValueError("Invalid mode. Expected one of: %s" % modes)

    # read in file
    # define intensity columns
    all_intcols, dict_numeric_cols = common.define_intensity_columns(sample1_colnames, sample2_colnames)

    # read in data
    df = common.read_data_table(file, dict_numeric_cols, all_intcols, pep_colname)

    # filter
    df_filt = common.filter_min_observed(df, grp1_intcols=sample1_colnames, grp2_intcols=sample2_colnames,
                                         threshold=threshold)

    results = 0
    descript = []
    if mode == 'fn':
        descript = ['name', 'namespace']
    if mode == 'taxfn':
        descript = ['descript']

    if mode == 'fn':
        results = functional_analysis(df=df_filt, go_colname=go_colname, all_intcols=all_intcols,
                                      grp1_intcols=sample1_colnames, grp2_intcols=sample2_colnames,
                                      test=test, threshold=threshold, outfile=outfile,
                                      ontology=ontology, slim_down=slim_down, paired=paired, obo_path=obo_path,
                                      slim_path=slim_path, download_obo=download_obo, overwrite_obo=overwrite_obo)


    if mode == 'tax':
        results = taxonomy_analysis(df=df_filt, all_intcols=all_intcols,
                                    sample1_colnames=sample1_colnames,
                                    sample2_colnames=sample2_colnames,
                                    test=test,
                                    paired=paired)

    if mode == 'taxfn':
        results = function_taxonomy_interaction_analysis(df=df_filt,
                                                         cog_name=cog_colname,
                                                         tax_rank=tax_rank,
                                                         all_intcols=all_intcols,
                                                         sample1_colnames=sample1_colnames,
                                                         sample2_colnames=sample2_colnames,
                                                         threshold=threshold,
                                                         outfile=outfile,
                                                         testtype=test_type,
                                                         paired=paired)

    if outfile:
        cols = []
        # order of columns we want
        if mode == 'tax':
            cols += ['member', 'rank']
        if test:
            cols += ['id', 'log2ratio_2over1', 'p', 'corrected_p']
            df_to_write = results.sort_values(by='corrected_p', axis=0, ascending=True)
        else:
            df_to_write = results
        if descript:
            cols += descript
        cols += all_intcols

        df_to_write.to_csv(outfile, sep="\t", header=True, index=False, columns=cols)

    return results


def main():
    pass

if __name__ == "__main__":
    main()