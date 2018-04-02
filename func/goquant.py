import goatools
import os
import pandas as pd


class ParentGoTerms:
    terms = {"biological_process":'GO:0008150',
             "molecular_function":'GO:0003674',
             "cellular_component":'GO:0005575'}
    namespaces = ['biological_process',
                  'molecular_function',
                  'cellular_component']


class GoQuant:
    def __init__(self, df, intcol):
        self.df = df
        self.go = goatools.obo_parser.GODag(os.getcwd() + '/data/go-basic.obo')
        self.intcol = intcol
        self.go_int = self.add_up_through_hierarchy()

    def set_of_all_parents(self, terms):
        all_rents = set(terms)
        for i in set(terms):
            all_rents.update(self.go[i].get_all_parents())
        return all_rents

    def add_up_through_hierarchy(self):
        go_dict = dict()
        # iterate through rows and assign intensity to all parents of each go term
        for index, row in self.df.iterrows():
            go_terms = row['go'].split(',')
            go_terms_parents = self.set_of_all_parents(go_terms)
            intensity = {x: row[x] for x in self.intcol}
            for term in go_terms_parents:
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
        namespace = [self.go[i].namespace for i in gos]

        # get dictionary with entry for each sample
        vals = {self.intcol[ints]: [go_dict[term][self.intcol[ints]] for term in gos] for ints in range(len(self.intcol))}
        namespace_dict = {"namespace": namespace}

        # merge the dictionaries
        merged_dict = {**vals, **namespace_dict}
        go_df = pd.DataFrame(merged_dict,
                             index=gos)
        # group by namespace, normalize to BP, MF, or CC
        namespaces = ParentGoTerms().namespaces
        for i in self.intcol:
            norm_name = "OSRA_" + i
            go_df[norm_name] = pd.concat([self.normalize_by_namespace(name, go_df, i) for name in namespaces])
        go_df['name'] = pd.Series([self.go.query_term(x).name for x in go_df.index], index=go_df.index)
        return go_df

    def normalize_by_namespace(self, namespace, go_df, intcol):
        parents = ParentGoTerms().terms
        if namespace in go_df['namespace'].values:
            go_parent_int = go_df.loc[parents[namespace]][intcol]
            subdf = go_df.loc[go_df.namespace == namespace, intcol]
            sub_norm=pd.Series(
                subdf.values / go_parent_int,
                index=subdf.index
            )
        else:
            sub_norm=pd.Series()
        return sub_norm

# TODO
# filereader
# write out
