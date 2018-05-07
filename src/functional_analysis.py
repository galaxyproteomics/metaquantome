import goatools
from goatools import mapslim
import pandas as pd
import wget
import os
import shutil
from src import common
import definitions

# TODO
# print unknown go terms
# support other functional terms (COG categories and E.C. numbers, for instance)


def functional_analysis(df, go_colname,
                        all_intcols,
                        grp1_intcols,
                        grp2_intcols,
                        test,
                        threshold,
                        outfile,
                        ontology,
                        slim_down,
                        paired,
                        obo_path,slim_path,
                        download_obo,
                        overwrite_obo):

    # read gos
    if not obo_path:
        obo_path = os.path.join(definitions.DATA_DIR, 'go', 'go-basic.obo')
    if not slim_path:
        slim_path = os.path.join(definitions.DATA_DIR, 'go', 'goslim_generic.obo')

    if download_obo:
        update_go_obo_file(obo_path,slim_path,overwrite_obo)

    go_dag = goatools.obo_parser.GODag(os.path.join(os.getcwd(), obo_path))

    if slim_down:
        go_dag_slim = goatools.obo_parser.GODag(os.path.join(os.getcwd(), slim_path))
    else:
        go_dag_slim = None

    # add up through hierarchy
    go_df = add_up_through_hierarchy(df, slim_down, go_dag, go_dag_slim, go_colname, all_intcols)

    # differential expression
    if test:
        # need to drop root terms to avoid NaN propagation
        go_df.drop(root_GO_terms.values(), inplace=True, errors="ignore")
        results = common.test_norm_intensity(go_df, grp1_intcols, grp2_intcols, paired)
    else:
        results = go_df

    return results


def add_up_through_hierarchy(df, slim_down, go_dag, go_dag_slim, gocol, all_intcols):

    # define the go ontology that is used to assign intensities to ancestors
    if slim_down:
        ref_go_dag = go_dag_slim
    else:
        ref_go_dag = go_dag

    go_dict = dict()
    # iterate through rows and assign intensity to all parents of each go term
    for index, row in df.iterrows():
        go_terms = row[gocol].split(',')
        go_terms_parents, unknown_gos = set_of_all_parents(go_terms, slim_down, go_dag, go_dag_slim)
        intensity = {x: row[x] for x in all_intcols}
        for term in go_terms_parents:
            if term in ref_go_dag.keys():
                if term in go_dict.keys():
                    current_term = go_dict[term]
                    for k, v in intensity.items():
                        current_term[k] += v
                else:
                    new_dict = dict()
                    for k, v in intensity.items():
                        new_dict[k] = v
                    go_dict[term] = new_dict

    # convert back to data frame
    gos = list(go_dict.keys())
    namespace = [ref_go_dag[i].namespace for i in gos]

    # get dictionary with entry for each sample
    vals = {all_intcols[ints]:
                [go_dict[term][all_intcols[ints]] for term in gos] for ints in range(len(all_intcols))}
    namespace_dict = {"namespace": namespace}

    # merge the dictionaries
    merged_dict = {**vals, **namespace_dict}
    go_df = pd.DataFrame(merged_dict,
                         index=gos)

    # group by namespace, normalize to BP, MF, or CC
    for i in all_intcols:
        norm_name = i
        go_df[norm_name] = pd.concat([normalize_by_namespace(name, go_df, i) for name in namespaces])
    go_df['name'] = pd.Series([ref_go_dag.query_term(x).name for x in go_df.index], index=go_df.index)
    return go_df


def set_of_all_parents(terms, slim_down, go_dag, go_dag_slim):
    unknown_gos = set()
    all_ancestors = set()
    if slim_down:
        for i in terms:
            if i in go_dag_slim.keys():
                # take all ancestors
                all_ancestors.update(mapslim.mapslim(i, go_dag, go_dag_slim)[1])
            else:
                unknown_gos.update([i])
    else:
        all_ancestors.update(terms)
        for i in set(terms):
            if i in go_dag.keys():
                all_ancestors.update(go_dag[i].get_all_parents())
            else:
                if i != "unknown":
                    unknown_gos.update([i])
    return all_ancestors, unknown_gos


def normalize_by_namespace(namespace, go_df, all_intcols):
    if namespace in go_df['namespace'].values:
        go_parent_int = go_df.loc[root_GO_terms[namespace]][all_intcols]
        sub_df = go_df.loc[go_df.namespace == namespace, all_intcols]
        sub_norm=pd.Series(
            sub_df.values / go_parent_int,
            index=sub_df.index
        )
    else:
        sub_norm=pd.Series()
    return sub_norm



namespaces = ['biological_process',
              'molecular_function',
              'cellular_component']
root_GO_terms = {"biological_process":'GO:0008150',
             "molecular_function":'GO:0003674',
             "cellular_component":'GO:0005575'}


def update_go_obo_file(obo_path, slim_path, overwrite=True):
    full_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    slim_obo_url = 'http://www.geneontology.org/ontology/subsets/goslim_generic.obo'
    full = wget.download(full_obo_url, out=obo_path)
    slim = wget.download(slim_obo_url, out=slim_path)
    if overwrite:
        if os.path.exists(obo_path):
            shutil.move(full, obo_path)
        if os.path.exists(slim_path):
            shutil.move(slim, slim_path)

