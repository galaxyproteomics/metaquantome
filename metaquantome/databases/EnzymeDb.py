import Bio.ExPASy.Enzyme as Enz
import os
import json
import re
import logging

from metaquantome.util.utils import stream_to_file_from_url


class EnzymeDb:
    # names to keep track of each level
    LEVEL_NAMES = ['ec0', 'ec1', 'ec2', 'ec3']

    # define all unknown as a list of 4 dashes.
    # this is consistent with a dash indicating that a
    # annotation is not known to that level of specificity
    ALL_UNKNOWN = ['-'] * 4

    def __init__(self, data_dir):
        """
        Create EnzymeDb object.

        :param data_dir: directory in which the ENZYME files will go
        """
        self.ecdb = self.load_enzyme_db(data_dir)

    @staticmethod
    def _define_ec_paths(data_dir):
        """
        define paths within data_dir corresponding to needed files
        for the ENZYME files

        :param data_dir: data directory
        :return: tuple of paths
        """
        dat_path = os.path.join(data_dir, 'enzyme.dat')
        dat_json = os.path.join(data_dir, 'ec_id.json')
        class_path = os.path.join(data_dir, 'enzclass.txt')
        class_json = os.path.join(data_dir, 'enzclass.json')

        return dat_path, dat_json, class_path, class_json

    @staticmethod
    def download_enzyme_db(data_dir, overwrite):
        """
        download enzyme db to files within data_dir

        :param data_dir: data directory
        :param overwrite: whether to overwrite existing files or not
        :return: None
        """
        dat_path, dat_json, class_path, class_json = EnzymeDb._define_ec_paths(data_dir)
        if (os.path.exists(dat_path) and os.path.exists(class_path)) and not overwrite:
            logging.info('Using ENZYME files in ' + data_dir)
        else:
            logging.info('Downloading enzyme files from ftp.expasy.org to ' + data_dir)
            enz_dat_url = 'ftp://ftp.expasy.org/databases/enzyme/enzyme.dat'
            stream_to_file_from_url(enz_dat_url, dat_path)
            enz_class_url = 'ftp://ftp.expasy.org/databases/enzyme/enzclass.txt'
            stream_to_file_from_url(enz_class_url, class_path)
            # create a simpler file from enzyme.dat
            EnzymeDb._create_ec_num_enzyme_name_association_file(dat_path, dat_json)
            # get clean enzclass.txt and dump to json
            EnzymeDb._read_enzyme_class_to_json(class_path, class_json)

    def load_enzyme_db(self, data_dir):
        """
        load enzyme database by loading the classes (3 or fewer levels)
        and the enzyme.dat file (4 levels of annotation)
        and merging them. This function calls _annotated_enzyme_db to
        annotate the terms.

        :param data_dir: data directory
        :return: annotated enzyme database
        """
        dat_path, dat_json, class_path, class_json = EnzymeDb._define_ec_paths(data_dir)
        with open(class_json, 'r') as enz_json:
            enz_class = json.load(enz_json)
        with open(dat_json) as ec:
            ec_id = json.load(ec)
        ecdb = {**enz_class, **ec_id}

        annotated_db = self._annotate_enzyme_db(ecdb)
        return annotated_db

    def _annotate_enzyme_db(self, ecdb):
        """
        Annotate the enzyme db with depth and description
        This function is split out mostly for testing.

        :param ecdb: non-annotated enzyme database
        :return: annotated enzyme database
        """
        annotated_db = dict()
        for k, v in ecdb.items():
            depth = self._assign_depth(k)
            levels = self._split_ec(k)
            descript = v
            newv = {'depth': depth,
                    'levels': levels,
                    'descript': descript,
                    'id': k}
            annotated_db[k] = newv
        return annotated_db

    @staticmethod
    def _create_ec_num_enzyme_name_association_file(enzdat_file, ec_id_file):
        """
        use the Enz class from BioPython to parse the enzyme.dat file,
        which contains all fully specified enzymes (i.e., with 4 levels)
        and annotations. The file is parsed to a dictionary, which is
        then dumped to JSON for easy reading later.

        :param enzdat_file: path to enzyme.dat file
        :param ec_id_file: path to JSON dump
        :return: None
        """
        with open(enzdat_file) as file:
            enzyme_dat = Enz.parse(file)
            enz_id = dict()
            for record in enzyme_dat:
                enz_id[record['ID']] = record['DE']
        with open(ec_id_file, 'w') as ec:
            json.dump(enz_id, ec)

    @staticmethod
    def _read_enzyme_class_to_json(class_path, class_json):
        """
        read the enzclass.txt file, which contains all levels of the ENZYME database higher than 4
        (i.e., it has 1.1.1.- but not 1.1.1.1).
        the file is parsed to get the class identifier and the description,
        which are inserted into a dictionary {class id: description, ...}

        the dictionary is then written to JSON, for quick reading later

        :param class_path: path to enzclass file
        :param class_json: path to write JSON to
        :return: None
        """
        enzyme_class = dict()
        with open(class_path) as f:
            for line in f:
                class_match = re.compile(r"^\d\.")
                id_pattern = re.compile(r"(^\d\.[ \a\d-]*\.[ \a\d-]*\.[ \a\d-]*)")
                descript_pattern = re.compile(r"([A-Z]{1}.*)$")
                if re.match(pattern=class_match, string=line[0:2]):
                    id = re.search(id_pattern, line)
                    clean_id = re.sub(pattern=r' ', repl='', string=id.group(1))
                    descript = re.search(descript_pattern, line).group(1)
                    enzyme_class[clean_id] = descript
        with open(class_json, 'w') as enz_json:
            json.dump(enzyme_class, enz_json)

    def is_in_db(self, ecid):
        """
        Check if the provided ecid is in the ENZYME database

        :param ecid: ecid as string
        :return: Boolean. if valid EC number, returns True.
        """
        if ecid in self.ecdb.keys():
            return True
        else:
            return False

    def get_children(self, ecid):
        """
        Get the children of the provided ec

        :param ecid: ecid as string
        :return: set of children. If no children, returns empty set.
        """
        parent = self.ecdb[ecid]
        parent_levels = parent['levels']
        parent_depth = parent['depth']

        # get annotation at levels up to depth
        annot_to_depth = parent_levels[0:(parent_depth+1)]

        # children are all terms that:
        #   - have same annotation at level 0 to parent_depth
        #   - are not the parent
        #   - have a depth that is one greater than the parent (i.e., immediate descendants)
        children = {ec for ec, annot in self.ecdb.items() if
                    annot['levels'][0:(parent_depth+1)] == annot_to_depth and
                    annot['depth'] == parent_depth + 1 and
                    ec != parent['id']}
        return children

    def get_descendants(self, ecid):
        """
        Get descendants of provided ecid

        :param ecid: ecid as string
        :return: set of all descendants of ecid. If no descendants, returns empty set
        """
        parent = self.ecdb[ecid]
        parent_levels = parent['levels']
        parent_depth = parent['depth']

        # get annotation at levels up to depth
        annot_to_depth = parent_levels[0:(parent_depth+1)]

        # descendants are all terms that:
        #   - have same annotation from level 0 to level parent_depth
        #   - are not the parent
        desc = {ec for ec, annot in self.ecdb.items() if
                annot['levels'][0:(parent_depth+1)] == annot_to_depth and
                ec != parent['id']}
        return desc

    def get_parents(self, ecid):
        """
        Returns the parents (immediate ancestors) of the provided ecid

        :param ecid: String of ecid in form 'w.x.y.z'
        :return: A set of all parents of the ecid. If no parents, returns empty set.
        """
        child = self.ecdb[ecid]
        child_levels = child['levels']
        child_depth = child['depth']
        # get annotation at levels up to depth - 1
        annot_to_prev_depth = child_levels[0:child_depth]
        # parents are all terms that:
        #   - have same annotation at level 0 to child_depth minus 1
        #   - have depth (child_depth - 1)
        #   - are not the child
        parents = {ec for ec, annot in self.ecdb.items() if
                   annot['levels'][0:child_depth] == annot_to_prev_depth and
                   annot['depth'] == child_depth - 1 and
                   ec != child['id']}
        return parents

    def get_ancestors(self, ecid):
        """
        Returns the ancestors of the provided ecid.

        :param ecid: String of ecid in form 'w.x.y.z'
        :return: A set of all ancestors of the ecid. If no ancestors, returns empty set.
        """
        # ancestors are all terms that:
        #   - have same annotation at child_depth
        #   - have depth less than child_depth
        #   - are not the child
        # run get_parents for each parent and each grandparent, etc.
        parents = self.get_parents(ecid)
        ancestors = parents.copy()
        while len(parents) > 0:
            this_id = parents.pop()
            this_parents = self.get_parents(this_id)
            ancestors.update(this_parents)
            parents.update(this_parents)
        return ancestors

    def _assign_depth(self, ecid):
        """
        Determine the depth of an EC number,
        where depth is defined as the most specific level indicated by the EC number.

        :param ecid: an enzyme classification number in the form 'w.x.y.z'.
        Unspecified levels must have a '-' (ex. 'w.x.-.-')
        :return: The depth of the ec number, ranging from 0 to 3.
        For example, 1.-.-.- has depth 0 and 1.2.3.4 has depth 3.
        """
        split_ec = self._split_ec(ecid)
        if '-' in split_ec:
            # get the first occurence of a '-', indicating unassigned at that level
            # subtract 1 to get the last assigned level
            depth = split_ec.index('-') - 1
        else:
            # if no '-' is present, that indicates maximum depth (3)
            depth = 3
        return depth

    @staticmethod
    def _split_ec(ecid):
        """
        split a string EC id into a list
        :param ecid: ecid as string
        :return: ecid as list
        """
        split_ec = ecid.split('.')
        return split_ec
