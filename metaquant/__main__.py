from metaquant import io,parser
from metaquant.runner import metaquant
import sys

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