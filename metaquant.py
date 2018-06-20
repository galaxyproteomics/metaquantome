from src.function_taxonomy_interaction import function_taxonomy_analysis
from src.functional_analysis import functional_analysis
from src.taxonomy_analysis import taxonomy_analysis
from src import io
from src import parser
import sys


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

def main():
    args = parser.parse_args_cli()

    # read sample names
    samp_grps = io.read_samp_info(args.samps)

    # input checks, error handling (todo)

    metaquant(args.mode, sample_names=samp_grps, int_file=args.int_file, pep_colname=args.pep_colname,
              func_file=args.func_file, tax_file=args.tax_file, ontology=args.ontology, tax_colname=args.tax_colname,
              outfile=args.outfile, slim_down=args.slim_down, test=args.test, paired=args.paired,
              threshold=args.threshold, obo_path=args.obo_path, slim_path=args.slim_path, update_obo=args.update_obo)

    sys.exit(0)


if __name__ == "__main__":
    main()