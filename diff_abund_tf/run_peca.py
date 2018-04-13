from rpy2.robjects.packages import importr
from rpy2.robjects.lib.dplyr import (DataFrame, mutate)
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import warnings

class TaxaFuncPECA:
    def __init__(self, file, cog_name, tax_name, list_samplename1, list_samplename2):
        mutate_string = "paste(" + cog_name + "," + tax_name + ", sep='-')"
        self.df = (DataFrame(DataFrame.from_csvfile(file, sep='\t', na_strings=["", "NA", "NaN"])).
              mutate(ft=mutate_string))
        self.samp1 = list_samplename1
        self.samp2 = list_samplename2

    def run_peca(self, test="t", paired=False, outfile=None):
        peca = importr("PECA")
        warnings.simplefilter('ignore')
        peca_results = peca.PECA_df(df=self.df, id="ft",
                                    samplenames1=ro.StrVector(self.samp1),
                                    samplenames2=ro.StrVector(self.samp2),
                                    test=test,
                                    paired=paired)
        peca_pandas = pandas2ri.ri2py_dataframe(peca_results)
        if outfile:
            peca_pandas.to_csv(outfile, sep='\t', index_label='function_taxon')
        return peca_pandas


# TODO
# add test type parameters, including test type and paired or not