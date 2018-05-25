import goatools
import os
from src import go
from src import common
import definitions
from src.cog import cogCat
from src.cog import take_first_cog

# TODO
# print unknown go terms
# add support for EC numbers


def functional_analysis(df, func_colname, all_intcols, grp1_intcols, grp2_intcols, test, threshold, ontology, slim_down,
                        paired, obo_path, slim_path, download_obo, overwrite_obo):

    if ontology == "go":
        # read gos
        if not obo_path:
            obo_path = os.path.join(definitions.DATA_DIR, 'go', 'go-basic.obo')
        if not slim_path:
            slim_path = os.path.join(definitions.DATA_DIR, 'go', 'goslim_generic.obo')

        if download_obo:
            go.update_go_obo_file(obo_path,slim_path,overwrite_obo)

        go_dag = goatools.obo_parser.GODag(os.path.join(os.getcwd(), obo_path))

        if slim_down:
            go_dag_slim = goatools.obo_parser.GODag(os.path.join(os.getcwd(), slim_path))
        else:
            go_dag_slim = None

        # add up through hierarchy
        go_df = go.add_up_through_hierarchy(df, slim_down, go_dag, go_dag_slim, func_colname, all_intcols)

        # differential expression
        if test:
            # need to drop root terms to avoid NaN propagation
            go_df.drop(go.ROOT_GO_TERMS.values(), inplace=True, errors="ignore")
            results = common.test_norm_intensity(go_df, grp1_intcols, grp2_intcols, threshold, paired)
        else:
            results = common.calc_means(go_df, grp1_intcols, grp2_intcols)

    elif ontology == "cog":
        cog_df = take_first_cog(df, func_colname)
        cog_sum_df = cog_df[[func_colname] + all_intcols].\
            groupby(func_colname).\
            sum()
        rel_cog_sum = cog_sum_df / cog_sum_df.sum(axis = 0)
        rel_cog_sum['descript'] = [cogCat[x] for x in rel_cog_sum.index]
        if test:
            results = common.test_norm_intensity(rel_cog_sum, grp1_intcols, grp2_intcols, threshold, paired, log=False)
        else:
            results = common.calc_means(cog_sum_df, grp1_intcols, grp2_intcols)

    else:
        raise ValueError("the desired ontology is not supported. " +
                         "Please use either GO (ontology = 'go') or COG (ontology = 'cog')")

    return results


