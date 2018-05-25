from src.function_taxonomy_interaction import function_taxonomy_interaction_analysis
from src.functional_analysis import functional_analysis
from src.taxonomy_analysis import taxonomy_analysis
from src import common

# TODO
# HIGH PRIORITY
# Allow an arbitrary number of samples (eventually, this could be in a text file)
# how to reconcile 2-sample tests with an arbitrary number of samples?
# one option is to state that testing can only be done with 2 samples, and only do if test is True and len(samps) = 2
# change to ncbi taxid input, because it is unambiguous
#   - translate to names at final step
#   - build dataframe and then do rank-specific normalization, as before
# LOWER PRIORITY
# write arg parser
# checks on input, error handling
#   - List of lists, where the inner lists are replicates and the outer lists indicate the experimental conditions
# Report group-wise mean (for all except TF)
# think more about TF

def metaquant(mode, file, sample1_colnames, sample2_colnames=None,
              cog_colname="cog",
              func_colname="go",
              lca_colname="lca",
              pep_colname="peptide",
              outfile=None, ontology="go",
              slim_down=False, test=False,
              test_type="modt",
              paired=False, threshold=0,
              obo_path=None, slim_path=None,
              download_obo=False, overwrite_obo=False):

    # read in file
    # define intensity columns
    all_intcols, dict_numeric_cols, grp1, grp2 = common.define_intensity_columns(sample1_colnames, sample2_colnames)

    # read in data
    df = common.read_data_table(file, dict_numeric_cols, all_intcols, pep_colname)

    # change ontology to lowercase
    lc_ontology = ontology.lower()

    modes = ['fn', 'tax', 'taxfn', 'fnpred']
    if mode == 'fn':
        results = functional_analysis(df=df, func_colname=func_colname, all_intcols=all_intcols,
                                      grp1_intcols=grp1, grp2_intcols=grp2, test=test,
                                      threshold=threshold, ontology=lc_ontology, slim_down=slim_down, paired=paired,
                                      obo_path=obo_path, slim_path=slim_path, download_obo=download_obo,
                                      overwrite_obo=overwrite_obo)

    elif mode == 'tax':
        results = taxonomy_analysis(df=df, all_intcols=all_intcols,
                                    sample1_colnames=grp1,
                                    sample2_colnames=grp2,
                                    test=test,
                                    threshold=threshold,
                                    paired=paired)

    elif mode == 'taxfn':
        results = function_taxonomy_interaction_analysis(df=df,
                                                         cog_name=cog_colname,
                                                         lca_colname=lca_colname,
                                                         sample1_colnames=grp1,
                                                         sample2_colnames=grp2,
                                                         threshold=threshold,
                                                         testtype=test_type,
                                                         paired=paired)
    else:
        raise ValueError("Invalid mode. Expected one of: %s" % modes)

    # descript column

    descript = []
    if mode == 'fn':
        if lc_ontology == "go":
            descript = ['name', 'namespace']
        if lc_ontology == "cog":
            descript = ['descript']

    if mode == 'taxfn':
        descript = ['descript']

    # outfile prettying
    if outfile:
        cols = ['id']
        if descript:
            cols += descript
        if mode != 'taxfn':
            if len(grp1) > 1:
                cols += ['mean1']
            if len(grp2) > 1:
                cols += ['mean2']

        # order of columns we want
        if mode == 'tax':
            cols += ['rank']
        if test:
            cols += ['log2ratio_2over1', 'p', 'corrected_p']
            df_to_write = results.sort_values(by='corrected_p', axis=0, ascending=True)
        else:
            df_to_write = results

        # in all cases but taxfn, we want the new intensity columns but we want them to be last
        if mode != 'taxfn':
            cols += all_intcols
        else:
            cols += ["n"]

        df_to_write.to_csv(outfile, sep="\t", header=True, index=False, columns=cols, na_rep="NA")

    return results


def main():
    pass

if __name__ == "__main__":
    main()