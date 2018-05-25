from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import warnings
from src.cog import cogCat
from src.cog import take_first_cog


def function_taxonomy_interaction_analysis(df, cog_name, lca_colname, samp_grps, testtype, paired):

    # take first cog
    df = take_first_cog(df, cog_name)

    pd_df = df.assign(ft=df[cog_name] + '-' + df[lca_colname])

    # keep only intensity columns and ft
    pd_df_int = pd_df[['ft'] + samp_grps.all_intcols]

    grp1 = samp_grps.grp_names[0]
    grp2 = samp_grps.grp_names[1]
    samples_in_grp1 = samp_grps.sample_names[grp1]
    samples_in_grp2 = samp_grps.sample_names[grp2]

    # transform data frame to R
    pandas2ri.activate()
    rdf = pandas2ri.py2ri(pd_df_int)

    # run peca
    peca = importr("PECA")
    warnings.simplefilter('ignore')
    # peca calculates ratios as 1 over 2
    peca_results = peca.PECA_df(df=rdf, id="ft",
                                samplenames1=ro.StrVector(samples_in_grp1),
                                samplenames2=ro.StrVector(samples_in_grp2),
                                test=testtype,
                                paired=paired)

    # translate back to pandas
    peca_pandas = pandas2ri.ri2py_dataframe(peca_results)
    ft = peca_pandas.index
    # add cog description
    peca_pandas['descript'] = [cogCat[x] for x in ft.str.split('-').str[0].values]

    # clearer column names, only keep important columns
    peca_pandas['id'] = peca_pandas.index

    fc_name = 'log2fc_' + grp1 + '_over_' + grp2
    peca_pandas.rename(index=str, columns={"slr": fc_name,
                                           "p.fdr": "corrected_p",
                                           "function_taxon": "id"}, inplace=True)
    peca_pandas.drop(columns=['t', 'score'], inplace=True)

    return peca_pandas
