import goatools
import os
import pandas as pd
import numpy as np


class ParentGoTerms:
    terms = {"biological_process":'GO:0008150',
             "molecular_function":'GO:0003674',
             "cellular_component":'GO:0005575'}
    namespaces = ['biological_process',
                  'molecular_function',
                  'cellular_component']

class GoQuant:
    def __init__(self, intcol, gocol, file):
        if not isinstance(intcol, list):
            self.intcol = [intcol]
        else:
            self.intcol = intcol
        self.gocol = gocol
        self.file = file
        self.df = self.read_table()
        self.go = goatools.obo_parser.GODag(os.getcwd() + '/data/go-basic.obo')
        self.unknown_gos = set()
        self.go_int = self.add_up_through_hierarchy()

        if self.unknown_gos == set():
            print("The following go terms were not found in the OBO file, and may be obsolete: ")
            for i in self.unknown_gos:
                print(i)

    def set_of_all_parents(self, terms):
        all_rents = set(terms)
        for i in set(terms):
            if i in self.go.keys():
                all_rents.update(self.go[i].get_all_parents())
            else:
                if i != "unknown":
                    self.unknown_gos.update([i])
        return all_rents

    def add_up_through_hierarchy(self):
        go_dict = dict()
        # iterate through rows and assign intensity to all parents of each go term
        for index, row in self.df.iterrows():
            go_terms = row[self.gocol].split(',')
            go_terms_parents = self.set_of_all_parents(go_terms)
            intensity = {x: row[x] for x in self.intcol}
            for term in go_terms_parents:
                if term in self.go.keys():
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

    def read_table(self):
        numeric_cols = {x: np.float64 for x in self.intcol}

        # read in data
        df = pd.read_table(self.file, sep="\t", index_col="peptide", dtype = numeric_cols)

        # drop columns where all are NA
        df.dropna(axis=1, how="all", inplace = True)

        # change missing intensities to 0
        df[self.intcol] = df[self.intcol].fillna(0)

        # change missing go annotation to 'unknown', because nan compares false to itself
        df.fillna('unknown', inplace=True)

        return df





# TODO
# write out
# add t.tests, fdr correction, etc. 
