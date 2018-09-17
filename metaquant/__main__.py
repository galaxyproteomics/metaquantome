from metaquant import runner
from metaquant.util import parser
import sys


def main():
    args = parser.parse_args_cli()
    runner.runner(args.mode, sinfo=args.samps, int_file=args.int_file, pep_colname=args.pep_colname,
                  func_file=args.func_file, tax_file=args.tax_file, ontology=args.ontology,
                  tax_colname=args.tax_colname, outfile=args.outfile, slim_down=args.slim_down, test=args.test,
                  paired=args.paired, threshold=args.threshold, data_dir=args.data_dir, overwrite=args.overwrite,
                  min_peptides=args.min_peptides, min_children_non_leaf=args.min_children_non_leaf)
    sys.exit(0)


if __name__ == "__main__":
    main()
