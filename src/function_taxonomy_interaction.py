from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import warnings


def function_taxonomy_interaction_analysis(df, cog_name, lca_colname,
                                           sample1_colnames, sample2_colnames, threshold,
                                           testtype, paired):

    # take first cog
    df[cog_name] = df[cog_name].str.split(',').str[0]

    pd_df = df.assign(ft=df[cog_name] + '-' + df[lca_colname])

    # keep only intensity columns and ft
    pd_df_int = pd_df[['ft'] + sample1_colnames + sample2_colnames]

    # # filter to only samples with certain threshold value
    # df_filt = common.filter_min_observed(pd_df_int, sample1_colnames, sample2_colnames, threshold)

    # transform data frame to R
    pandas2ri.activate()
    rdf = pandas2ri.py2ri(pd_df_int)

    # run peca
    peca = importr("PECA")
    warnings.simplefilter('ignore')
    # peca calculates ratios as 1 over 2, so we exchange the samples
    peca_results = peca.PECA_df(df=rdf, id="ft",
                                samplenames2=ro.StrVector(sample1_colnames),
                                samplenames1=ro.StrVector(sample2_colnames),
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

# from ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab

cogCat = {'J': 'Translation, ribosomal structure and biogenesis',
          'A': 'RNA processing and modification',
          'K': 'Transcription',
          'L': 'Replication, recombination and repair',
          'B': 'Chromatin structure and dynamics',
          'D': 'Cell cycle control, cell division, chromosome partitioning',
          'Y': 'Nuclear structure',
          'V': 'Defense mechanisms',
          'T': 'Signal transduction mechanisms',
          'M': 'Cell wall/membrane/envelope biogenesis',
          'N': 'Cell motility',
          'Z': 'Cytoskeleton',
          'W': 'Extracellular structures',
          'U': 'Intracellular trafficking, secretion, and vesicular transport',
          'O': 'Posttranslational modification, protein turnover, chaperones',
          'X': 'Mobilome: prophages, transposons',
          'C': 'Energy production and conversion',
          'G': 'Carbohydrate transport and metabolism',
          'E': 'Amino acid transport and metabolism',
          'F': 'Nucleotide transport and metabolism',
          'H': 'Coenzyme transport and metabolism',
          'I': 'Lipid transport and metabolism',
          'P': 'Inorganic ion transport and metabolism',
          'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
          'R': 'General function prediction only',
          'S': 'Function unknown',
          'unknown': 'unassigned function'}
