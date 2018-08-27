from goatools import mapslim
from goatools import obo_parser
import pandas as pd
import wget
import os
from metaquant.utils import safe_cast_to_list
import logging


class GeneOntologyDb:
    NAMESPACES = ['biological_process',
                  'molecular_function',
                  'cellular_component']
    ROOT_GO_TERMS = {"biological_process":'GO:0008150',
                     "molecular_function":'GO:0003674',
                     "cellular_component":'GO:0005575'}

    def __init__(self, data_dir, slim_down=False, overwrite=False):
        gofull, goslim = self._go_database_handler(data_dir, slim_down, overwrite)
        self.gofull = gofull
        self.goslim = goslim

    @staticmethod
    def _go_database_handler(data_dir, slim_down, overwrite):
        obo_path = os.path.join(data_dir, 'go-basic.obo')
        slim_path = os.path.join(data_dir, 'goslim_generic.obo')
        if (os.path.exists(obo_path) or os.path.exists(slim_path)) and not overwrite:
            logging.info('Using GO files in ' + data_dir)
        else:
            full_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
            logging.info('Downloading full GO obo file from ' + full_obo_url + ' to ' + obo_path)
            wget.download(full_obo_url, out=obo_path)
            slim_obo_url = 'http://www.geneontology.org/ontology/subsets/goslim_generic.obo'
            logging.info('Downloading generic slim GO obo file from ' + slim_obo_url + ' to ' + slim_path)
            wget.download(slim_obo_url, out=slim_path)
        # read gos
        go_dag = obo_parser.GODag(obo_path)
        if slim_down:
            go_dag_slim = obo_parser.GODag(slim_path)
        else:
            go_dag_slim = None
        return go_dag, go_dag_slim

    def map_to_slim(self):
        # todo figure out how to implement this with get_* as written
        # maybe want the ancestor that is closest to the original term by path?
        pass

    def _safe_query_go(self, goid):
        if goid in self.gofull.keys():
            return True, self.gofull[goid]
        else:
            return False, set()

    def get_children(self, goid):
        is_in_go, term = self._safe_query_go(goid)
        if is_in_go:
            children = {child.id for child in term.children}
            return children
        else:
            return set()

    def get_descendants(self, goid):
        is_in_go, term = self._safe_query_go(goid)
        if is_in_go:
            children = set(term.children)
            descendants = children.copy()
            while len(children) > 0:
                this_term = children.pop()
                this_children = this_term.children
                descendants.update(this_children)
                children.update(this_children)
            desc_ids = {term.id for term in descendants}
            return desc_ids
        else:
            return set()

    def get_parents(self, goid):
        is_in_go, term = self._safe_query_go(goid)
        if is_in_go:
            parents = {parent.id for parent in term.parents}
            return parents
        else:
            return set()

    def get_ancestors(self, goid):
        is_in_go, term = self._safe_query_go(goid)
        if is_in_go:
            parents = set(term.parents)
            ancestors = parents.copy()
            while len(parents) > 0:
                this_term = parents.pop()
                this_parents = this_term.parents
                ancestors.update(this_parents)
                parents.update(this_parents)
            anc_ids = {term.id for term in ancestors}
            return anc_ids
        else:
            return set()


#
# def add_up_through_hierarchy(df, slim_down, go_dag, go_dag_slim, gocol, all_intcols, norm_to_root=False):
#     # drop peptides without go annotation
#     df.dropna(axis=0, how='any', subset=safe_cast_to_list(gocol), inplace=True)
#
#     # define the go ontology that is used to assign intensities to ancestors
#     if slim_down:
#         ref_go_dag = go_dag_slim
#     else:
#         ref_go_dag = go_dag
#
#     go_dict = dict()
#
#     # iterate through rows and assign intensity to all parents of each go term
#     # todo: keep track of number of peptides here
#     for index, row in df.iterrows():
#         go_terms = row[gocol].split(',')
#         go_terms_parents, unknown_gos = set_of_all_parents(go_terms, slim_down, go_dag, go_dag_slim)
#         intensity = {x: row[x] for x in all_intcols}
#         for term in go_terms_parents:
#             if term in ref_go_dag.keys():
#                 if term in go_dict.keys():
#                     current_term = go_dict[term]
#                     for k, v in intensity.items():
#                         current_term[k] += v
#                 else:
#                     new_dict = dict()
#                     for k, v in intensity.items():
#                         new_dict[k] = v
#                     go_dict[term] = new_dict
#
#     # convert back to data frame
#     gos = list(go_dict.keys())
#     namespace = [ref_go_dag[i].namespace for i in gos]
#
#     # get dictionary with entry for each sample
#     vals = {all_intcols[ints]:
#                 [go_dict[term][all_intcols[ints]] for term in gos] for ints in range(len(all_intcols))}
#     namespace_dict = {"namespace": namespace}
#
#     # merge the dictionaries
#     merged_dict = {**vals, **namespace_dict}
#     go_df = pd.DataFrame(merged_dict,
#                          index=gos)
#
#     # group by namespace, normalize to BP, MF, or CC
#     if norm_to_root:
#         for i in all_intcols:
#             norm_name = i
#             go_df[norm_name] = pd.concat([normalize_by_namespace(name, go_df, i) for name in self.ROOT_GO_TERMS.keys()])
#     go_df['id'] = go_df.index
#     go_df['name'] = pd.Series([ref_go_dag.query_term(x).name for x in go_df.index], index=go_df.index)
#     return go_df
#
#
# def set_of_all_parents(terms, slim_down, go_dag, go_dag_slim):
#     unknown_gos = set()
#     all_ancestors = set()
#     if slim_down:
#         for i in terms:
#             if i in go_dag.keys():
#                 # take all ancestors
#                 all_ancestors.update(mapslim.mapslim(i, go_dag, go_dag_slim)[1])
#             else:
#                 unknown_gos.update([i])
#     else:
#         all_ancestors.update(terms)
#         for i in set(terms):
#             if i in go_dag.keys():
#                 all_ancestors.update(go_dag[i].get_all_parents())
#             else:
#                 if i != "unknown":
#                     unknown_gos.update([i])
#     return all_ancestors, unknown_gos
#
#
# def normalize_by_namespace(namespace, go_df, all_intcols):
#     if namespace in go_df['namespace'].values:
#         go_parent_int = go_df.loc[ROOT_GO_TERMS[namespace]][all_intcols]
#         sub_df = go_df.loc[go_df.namespace == namespace, all_intcols]
#         sub_norm=pd.Series(
#             sub_df.values / go_parent_int,
#             index=sub_df.index
#         )
#     else:
#         sub_norm=pd.Series()
#     return sub_norm
#


