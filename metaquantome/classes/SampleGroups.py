import itertools
import json
import os
import numpy as np


class SampleGroups:
    """
    A class that contains information about which columns in the intensity data frame correspond to
    which experimental conditions.
    Other members of the class include the output column names for intensity, number of sample children, and
    number of unique peptides.
    """
    def __init__(self, sinfo):
        """
        Build SampleGroups object using sample info

        :param sinfo: Either a JSON-formatted string or a file path - see read_samp_info
        """
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

        # this is used in the stat module
        self.fc_name = None
        if self.ngrps == 2:
            grp1 = self.grp_names[0]
            grp2 = self.grp_names[1]
            self.fc_name = 'log2fc_' + grp1 + '_over_' + grp2

    def read_samp_info(self, sinfo):
        """
        read sample object to a dictionary

        :param sinfo: either a path to a tabular file or a JSON-formatted string
        :return: a dictionary of the form {group1_name: [grp1_samp1, grp1_samp2, ...], ...}
        """
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

    @staticmethod
    def to_json(obj):
        """
        checks that an object is correctly JSON formatted, and parses it if so
        thanks to https://stackoverflow.com/questions/5508509/how-do-i-check-if-a-string-is-valid-json-in-python

        :param obj: string
        :return: tuple. First element is True if JSON, False if not.
        Second element is json object (if JSON) or the original object (if not)
        """
        try:
            json_object = json.loads(obj)
            return True, json_object
        except ValueError:
            return False, obj
