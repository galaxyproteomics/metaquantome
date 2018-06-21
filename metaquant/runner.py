from metaquant.function_taxonomy_interaction import function_taxonomy_analysis
from metaquant.functional_analysis import functional_analysis
from metaquant.taxonomy_analysis import taxonomy_analysis
from metaquant import io


def metaquant(mode, sample_names, int_file, pep_colname='peptide', func_file=None, tax_file=None, ontology='go',
              tax_colname=None, outfile=None, slim_down=False, test=False, paired=False, threshold=0, obo_path=None,
              slim_path=None, update_obo=False):

    # read in file
    # define object with sample groups, intensity columns, etc.
    samp_grps = io.SampleGroups(sample_names)

    modes = ['fn', 'tax', 'taxfn']
    if mode == 'fn':
        df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                    int_file=int_file, func_file=func_file, func_colname=ontology)
        results = functional_analysis(df=df, func_colname=ontology, samp_grps=samp_grps, test=test, threshold=threshold,
                                      ontology=ontology, slim_down=slim_down, paired=paired, obo_path=obo_path,
                                      slim_path=slim_path, download_obo=update_obo)

    elif mode == 'tax':
        df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                    int_file=int_file, tax_file=tax_file, tax_colname=tax_colname)
        results = taxonomy_analysis(df=df, samp_grps=samp_grps,
                                    test=test,
                                    threshold=threshold,
                                    paired=paired,
                                    tax_colname=tax_colname)

    elif mode == 'taxfn':
        if ontology == 'cog':
            cog_colname = 'cog'
        else:
            raise ValueError("Only cog is supported for ft interaction. Make sure you have a cog column, named 'cog'")

        df = io.read_and_join_files(mode, pep_colname, samp_grps,
                                    int_file=int_file,
                                    tax_file=tax_file, tax_colname=tax_colname,
                                    func_file=func_file, func_colname=cog_colname)
        results = function_taxonomy_analysis(df=df, cog_name=cog_colname, lca_colname=tax_colname,
                                             samp_grps=samp_grps, test=test, threshold=threshold,
                                             paired=paired)

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


