from goatools import obo_parser
import os
import logging

from metaquantome.util.utils import safe_cast_to_list, stream_to_file_from_url


class GeneOntologyDb:
    # the three "namespaces", or "ontologies"
    NAMESPACES = ['biological_process',
                  'molecular_function',
                  'cellular_component']

    # the root node (most general) of each ontology
    ROOT_GO_TERMS = {"biological_process": 'GO:0008150',
                     "molecular_function": 'GO:0003674',
                     "cellular_component": 'GO:0005575'}

    def __init__(self, data_dir, slim_down=False):
        """
        create GeneOntologyDb object

        :param data_dir: directory that contains the database files
        :param slim_down: whether slim will be used in the analysis. If False, the slim db is not loaded
        (self.goslim = None)
        """
        gofull, goslim = self._load_go_db(data_dir, slim_down)
        self.gofull = gofull
        self.goslim = goslim
        self.slim_down = slim_down
        # slim_members is a set of all terms in the GO slim
        self.slim_members = None
        if slim_down:
            self.slim_members = set(goslim.keys())

    @staticmethod
    def _define_data_paths(data_dir):
        """
        define OBO and slim OBO paths within data dir.

        :param data_dir: Data directory
        :return: path to the full OBO and the slim OBO, as a tuple
        """
        obo_path = os.path.join(data_dir, 'go-basic.obo')
        slim_path = os.path.join(data_dir, 'goslim_metagenomics.obo')
        return obo_path, slim_path

    @staticmethod
    def download_go(data_dir, overwrite):
        """
        Download both full GO and GO metagenomics slim to data_dir.

        :param data_dir: Data directory
        :param overwrite: Whether to overwrite existing databases or not.
        :return: None
        """
        obo_path, slim_path = GeneOntologyDb._define_data_paths(data_dir)
        if (os.path.exists(obo_path) and os.path.exists(slim_path)) and not overwrite:
            logging.info('GO files exist in specified directory ' +
                         'and --update was not provided. Doing nothing.')
        else:
            full_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
            logging.info('Downloading full GO obo file from ' + full_obo_url + ' to ' + obo_path)
            stream_to_file_from_url(full_obo_url, obo_path)

            slim_obo_url = 'http://current.geneontology.org/ontology/subsets/goslim_metagenomics.obo'
            logging.info('Downloading generic slim GO obo file from ' + slim_obo_url + ' to ' + slim_path)
            stream_to_file_from_url(slim_obo_url, slim_path)

    @staticmethod
    def _load_go_db(data_dir, slim_down):
        """
        Load GO databases using goatools.oboparser. Always
        loads the full GO, and loads the metagenomics slim GO if slim_down = True

        :param data_dir: Data directory
        :param slim_down: Whether slim database is going to be used or not
        :return: A tuple, with full and slim GO. If slim_down = False, then the tuple is (full GO, None)
        """
        obo_path, slim_path = GeneOntologyDb._define_data_paths(data_dir)
        if not (os.path.exists(obo_path) and os.path.exists(slim_path)):
            logging.error('GO files not found in specified directory.\n' +
                          'Please use the command >metaquantome db ...  to download the files.')
        # read gos
        go_dag = obo_parser.GODag(obo_path)
        if slim_down:
            go_dag_slim = obo_parser.GODag(slim_path)
        else:
            go_dag_slim = None
        return go_dag, go_dag_slim

    def is_in_db(self, goid):
        """
        test whether a specific term is in the full GO database

        :param goid: query GO term
        :return: True if GO is present in GO database; False otherwise
        """
        if goid in self.gofull.keys():
            return True
        else:
            return False

    def map_set_to_slim(self, sample_set):
        """
        Maps a set of GO terms to the generic GO slim.
        For each term, the closest ancestor contained within the GO slim is reported.
        The return type is a dict, to account for multiple terms mapping to the same slim.
        See map_id_to_slim() for more details.

        :param sample_set: set of all GO terms reported by the annotation tool
        :return: dictionary where the keys are the members of sample_set and the values are the mapped-to slim terms
        """
        mapper = dict()
        for id in sample_set:
            mapper[id] = self.map_id_to_slim(id)
        return mapper

    def map_id_to_slim(self, goid):
        """
        Maps a GO term ID to the most closely related term in the generic GO slim.
        If goid is in GO slim, it is returned. If it is not even in the full GO, the string "unknown" is returned.
        Otherwise, the most recent ancestor in the GO slim is returned.
        If multiple ancestors tie, the one that is first alphabetically is chosen.

        :param goid: GO term ID
        :return: closest term in slim (which is the term itself if it is in the slim GO)
        """
        # first, check that term is in full GO. if not, return "unknown"
        if not self._safe_query_go(goid):
            return "unknown"
        # first, see if term is in slim
        slim_ids = self.goslim.keys()
        if goid in slim_ids:
            return goid
        # if not in slim, get closest ancestor that is in slim. if there is a tie,
        # select term that is first alphabetically (todo: could change this)
        # this should always finish because each of BP, CC, and MF are in slim
        closest_in_slim = None
        ancestor_set = self.get_parents(goid)
        while not closest_in_slim:
            potential_closest = set()
            # get term parents
            ancestors = list(ancestor_set)
            in_slim = [full_id in slim_ids for full_id in ancestors]
            for i in range(0, len(ancestors)):
                if in_slim[i]:
                    potential_closest.update(safe_cast_to_list(ancestors[i]))
            # no potential closest
            if len(potential_closest) == 0:
                # new ancestors are the parents of the old ancestors
                parents_of_new_ancestors = set()
                for term in ancestors:
                    parents_of_new_ancestors.update(self.get_parents(term))
                ancestor_set = parents_of_new_ancestors
                continue
            # one or more potential closest
            else:
                potential_ancestors_list = list(potential_closest)
                first_alpha = sorted(potential_ancestors_list)[0]
                closest_in_slim = first_alpha
        return closest_in_slim

    def _safe_query_go(self, goid):
        """
        normally, a dict throws an error if the requested key does not exist.
        This overrides that behavior so that None is returned

        :param goid: query GO id
        :return: the GO id object from the GO database, or None
        """
        if goid in self.gofull.keys():
            return self.gofull[goid]
        else:
            return None

    def get_children(self, goid):
        """
        get the immediate children of a go term

        :param goid: query GO id
        :return: children of the GO term. An empty set is returned
        if the term is not present in the database
        or it has no children.
        """
        term = self._safe_query_go(goid)
        if term:
            children = {child.id for child in term.children}
            return children
        else:
            return set()

    def get_descendants(self, goid):
        """
        Get all descendants of a GO term.

        :param goid: query GO id
        :return: set of descendants or empty set, if term doesn't exist in db or has no descendants
        """
        term = self._safe_query_go(goid)
        if term:
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
        """
        Get parents (immediate ancestors) of a GO term.

        :param goid: query GO id
        :return: set of parents or empty set, if term doesn't exist in db or has no parents
        """
        term = self._safe_query_go(goid)
        if term:
            parents = {parent.id for parent in term.parents}
            return parents
        else:
            return set()

    def get_ancestors(self, goid):
        """
        Get all ancestors of a GO term

        :param goid: query GO id
        :return: set of ancestors or empty set, if term doesn't exist in db or has no ancestors
        """
        term = self._safe_query_go(goid)
        if term:
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
