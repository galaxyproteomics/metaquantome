from src import go
from src import stats
from src.cog import cogCat
from src.cog import take_first_cog
import numpy as np


def functional_analysis(df, func_colname, samp_grps, test, threshold, ontology, slim_down, paired, obo_path, slim_path,
                        download_obo):

    if ontology == "go":
        go_dag, go_dag_slim = go.load_obos(obo_path, slim_path, slim_down, download_obo)

        # add up through hierarchy
        df_to_return = go.add_up_through_hierarchy(df, slim_down,
                                                   go_dag, go_dag_slim, func_colname, samp_grps.all_intcols)
        df_to_return.drop(go.ROOT_GO_TERMS.values(), inplace=True, errors="ignore")

    elif ontology == "cog":
        cog_df = take_first_cog(df, func_colname)
        cog_sum_df = cog_df[[func_colname] + samp_grps.all_intcols].\
            groupby(func_colname).\
            sum()

        df_to_return = cog_sum_df

        df_to_return['description'] = [cogCat[x] for x in df_to_return.index]

    else:
        raise ValueError("the desired ontology is not supported. " +
                         "Please use either GO (ontology = 'go') or COG (ontology = 'cog')")

    # differential expression
    if test and samp_grps.ngrps == 2:
        results = stats.test_norm_intensity(df_to_return, samp_grps, threshold, paired, log=False)
    else:
        results = stats.calc_means(df_to_return, samp_grps)

    # take log of intensities for return
    # replace zeros with nan
    results[results == 0] = np.nan

    results[samp_grps.all_intcols] = np.log2(results[samp_grps.all_intcols])

    # go term id numbers
    if ontology == 'go':
        results['go_id'] = results.index
    if ontology == 'cog':
        results['cog'] = results.index

    return results


