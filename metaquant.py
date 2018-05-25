from src.function_taxonomy_interaction import function_taxonomy_interaction_analysis
from src.functional_analysis import functional_analysis
from src.taxonomy_analysis import taxonomy_analysis
from src import common

def metaquant(mode, file, sample_names,
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
    # define object with sample groups, intensity columns, etc.
    samp_grps = common.SampleGroups(sample_names)

    # read in data
    df = common.read_data_table(file, samp_grps, pep_colname)

    # change ontology to lowercase
    lc_ontology = ontology.lower()

    modes = ['fn', 'tax', 'taxfn', 'fnpred']
    if mode == 'fn':
        results = functional_analysis(df=df, func_colname=func_colname, samp_grps=samp_grps, test=test,
                                      threshold=threshold, ontology=lc_ontology, slim_down=slim_down, paired=paired,
                                      obo_path=obo_path, slim_path=slim_path, download_obo=download_obo,
                                      overwrite_obo=overwrite_obo)

    elif mode == 'tax':
        results = taxonomy_analysis(df=df, samp_grps=samp_grps, test=test, threshold=threshold, paired=paired)

    elif mode == 'taxfn':
        results = function_taxonomy_interaction_analysis(df=df, cog_name=cog_colname, lca_colname=lca_colname,
                                                         samp_grps=samp_grps, testtype=test_type, paired=paired)
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
            cols += samp_grps.mean_names

        # order of columns we want
        if mode == 'tax':
            cols += ['rank']
        if test:
            fc_name = 'log2fc_' + samp_grps.grp_names[0] + '_over_' + samp_grps.grp_names[1]
            cols += [fc_name, 'p', 'corrected_p']
            df_to_write = results.sort_values(by='corrected_p', axis=0, ascending=True)
        else:
            df_to_write = results

        # in all cases but taxfn, we want the new intensity columns but we want them to be last
        if mode != 'taxfn':
            cols += samp_grps.all_intcols

        else:
            cols += ["n"]

        df_to_write.to_csv(outfile, sep="\t", header=True, index=False, columns=cols, na_rep="NA")

    return results


def main():
    pass

if __name__ == "__main__":
    main()