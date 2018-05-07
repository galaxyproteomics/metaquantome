from rpy2.robjects.packages import importr
from rpy2.robjects.lib.dplyr import (DataFrame, mutate)
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import warnings
import pandas as pd
from src import common

# TODO
# what to do about missing cog and missing taxon


def function_taxonomy_interaction_analysis(df, cog_name, tax_rank, all_intcols,
                                           sample1_colnames, sample2_colnames,
                                           outfile, threshold,
                                           testtype, paired):

    # take first cog
    df[cog_name] = df[cog_name].str.split(',').str[0]

    pd_df = df.assign(ft=df[cog_name] + '-' + df[tax_rank])

    # transform data frame to R
    pandas2ri.activate()
    rdf = pandas2ri.py2ri(pd_df)

    # run peca
    peca = importr("PECA")
    warnings.simplefilter('ignore')
    peca_results = peca.PECA_df(df=rdf, id="ft",
                                samplenames1=ro.StrVector(sample1_colnames),
                                samplenames2=ro.StrVector(sample2_colnames),
                                test=testtype,
                                paired=paired)

    # translate back to pandas
    peca_pandas = pandas2ri.ri2py_dataframe(peca_results)
    ft = peca_pandas.index
    # add cog description
    peca_pandas['descript'] = [cogCat[x] for x in ft.str.split('-').str[0].values]

    # clearer column names, only keep important columns
    peca_pandas['id'] = peca_pandas.index
    peca_pandas.rename(index=str, columns={"slr": "log2ratio_2over1",
                                           "p.fdr": "corrected_p",
                                           "function_taxon": "id"}, inplace=True)
    peca_pandas.drop(columns=['t', 'score'], inplace=True)

    return peca_pandas


cogCat = {'A': 'RNA processing and modification',
          'B': 'Chromatin structure and dynamics',
          'C': 'Energy production and conversion',
          'D': 'Cell cycle control and mitosis',
          'E': 'Amino acid metabolism and transport',
          'F': 'Nucleotide metabolism and transport',
          'G': 'Carbohydrate metabolism and transport',
          'H': 'Coenzyme metabolism',
          'I': 'Lipid metabolism',
          'J': 'Translation',
          'K': 'Transcription',
          'L': 'Replication and repair',
          'M': 'Cell wall/membrane/envelope biogenesis',
          'N': 'Cell motility',
          'O': 'Post-translational modification, protein turnover, chaperone functions',
          'P': 'Inorganic ion transport and metabolism',
          'Q': 'Secondary structure',
          'T': 'Signal transduction',
          'U': 'Intracellular trafficking and secretion',
          'Y': 'Nuclear structure',
          'Z': 'Cytoskeleton',
          'R': 'General Function Prediction only',
          'S': 'Function unknown',
          'nan': 'Function + taxa unknown'}


