import itertools
import json
import os
import numpy as np


class SampleGroups:
    def __init__(self, sinfo):
        # todo: doc
        # top level dictionary
        sample_names = self.read_samp_info(sinfo)
        self.sample_names = sample_names

        # flatten sample column names
        self.all_intcols = sorted(list(itertools.chain(*list(sample_names.values()))))

        # for defining pandas column data types on read in
        self.dict_numeric_cols = {x: np.float64 for x in self.all_intcols}

        # number of conditions
        self.ngrps = len(sample_names)

        # name of experimental groups
        # sort alphabetically, so it's deterministic
        self.grp_names = sorted(list(sample_names.keys()))

        # when calculating means, column names for means
        # same order as grp names
        self.mean_names = [grp + "_mean" for grp in self.grp_names]

        # other quant names
        samp_child_name = '_n_samp_children'
        self.samp_children_names_dict = {grp: [samp + samp_child_name for samp in self.sample_names[grp]]
                                         for grp in self.sample_names.keys()}
        self.samp_children_names_flat = [samp + samp_child_name for samp in self.all_intcols]

        samp_npep_name = '_n_peptide'
        self.n_peptide_names_dict = {grp: [samp + samp_npep_name for samp in self.sample_names[grp]]
                                     for grp in self.sample_names.keys()}
        self.n_peptide_names_flat = [samp + samp_npep_name for samp in self.all_intcols]

        # all columns that should be numeric when reading expanded
        expand_num = self.all_intcols +\
                     self.samp_children_names_flat +\
                     self.n_peptide_names_flat
        self.dict_numeric_cols_expanded = {x: np.float64 for x in expand_num}

        # this is used in test()
        self.fc_name = None
        if self.ngrps == 2:
            grp1 = self.grp_names[0]
            grp2 = self.grp_names[1]
            self.fc_name = 'log2fc_' + grp1 + '_over_' + grp2

    def read_samp_info(self, sinfo):
        # todo: doc
        # check if sinfo is a file name
        if os.path.exists(sinfo):
            samp_names = dict()
            with open(sinfo, 'r') as f:
                # throw away header
                f.readline()

                # read groups one by one
                for line in f:
                    split = line.split('\t')
                    samp_names[split[0]] = [elem.strip() for elem in split[1].split(',')]
            return samp_names
        else:
            # check if sinfo is a json format
            is_json, json_obj = self.to_json(sinfo)
            if is_json:
                return json_obj
            else:
                raise ValueError('--samps is not a text file or in proper json format. please check again!')

    # thanks to https://stackoverflow.com/questions/5508509/how-do-i-check-if-a-string-is-valid-json-in-python
    def to_json(self, obj):
        # todo: doc
        try:
            json_object = json.loads(obj)
            return True, json_object
        except ValueError:
            return False, obj