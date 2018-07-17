from metaquant import go
from metaquant import stats
from metaquant.cog import cogCat
from metaquant.cog import take_first_cog
import metaquant.ec as ec
import numpy as np
import pandas as pd
from metaquant.definitions import DATA_DIR
import os


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

    elif ontology == "ec":
        ec.enzyme_database_handler(download_obo, os.path.join(DATA_DIR, 'enzyme'))

        ec_df = df['ec'].apply(ec.expand_ec).join(df)

        # new - just add up through ranks
        ec_sum_df = pd.concat(
            [ec.abundance_level(ec_df, x, samp_grps.all_intcols, norm_to_rank=False) for x in ec.LEVEL_NAMES])

        df_to_return = ec_sum_df

        ec_descript_dict = ec.load_combined_enzyme_class_ec_id(os.path.join(DATA_DIR, 'enzyme'))
        df_to_return['description'] =\
            [ec_descript_dict[x] if x in ec_descript_dict.keys() else 'unknown_ec' for x in df_to_return.index]

    else:
        raise ValueError("the desired ontology is not supported. " +
                         "Please use either GO (ontology = 'go'), COG (ontology = 'cog'), or EC numbers (ontology = 'ec')")

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


