from ete3 import NCBITaxa
import os
import numpy as np
import logging
import warnings


# each of the databases should have these methods:
# 1. get_children
# 2. get_parents
# 3. get_ancestors
# 4. get_descendants
# more?

BASIC_TAXONOMY_TREE = ["phylum",
                       "class",
                       "order",
                       "family",
                       "genus",
                       "species"]
FULL_TAXONOMY_TREE = ['no rank',
                      'superkingdom',
                      'kingdom',
                      'subkingdom',
                      'superphylum',
                      'phylum',
                      'subphylum',
                      'superclass',
                      'class',
                      'subclass',
                      'infraclass',
                      'superorder',
                      'order',
                      'suborder',
                      'infraorder',
                      'parvorder',
                      'superfamily',
                      'family',
                      'subfamily',
                      'tribe',
                      'subtribe',
                      'genus',
                      'subgenus',
                      'species group',
                      'species subgroup',
                      'species',
                      'subspecies',
                      'varietas',
                      'forma']

NUMERIC_RANK = {FULL_TAXONOMY_TREE[i]: i for i in range(len(FULL_TAXONOMY_TREE))}
BASIC_NUMERIC_RANK = [NUMERIC_RANK[i] for i in BASIC_TAXONOMY_TREE]
UNKNOWN_RANK = 'unknown'
UNIDENTIFIED = 32644


class NCBITaxonomyDb:
    def __init__(self, data_dir):
        # todo: doc
        self.ncbi = self._load_ncbi_db(data_dir)

    @staticmethod
    def _define_tax_paths(data_dir):
        return os.path.join(data_dir, 'taxa.sqlite')

    @staticmethod
    def download_ncbi(data_dir):
        tax_path = NCBITaxonomyDb._define_tax_paths(data_dir)
        if os.path.exists(tax_path):
            logging.info('Database exists. Skipping download. Database is automatically updated each time used.')
        else:
            with warnings.catch_warnings():  # turning off ResourceWarnings (unclosed file) from ete3
                warnings.simplefilter('ignore')
                open(tax_path, 'a').close()  # will say that database is not up to date
                NCBITaxa(tax_path)  # this prints a lot of logging messages, so no message here
                # remove taxdump, in working directory
                taxdump = os.path.join(os.getcwd(), 'taxdump.tar.gz')
                if os.path.exists(taxdump):
                    os.remove(taxdump)

    @staticmethod
    def _load_ncbi_db(data_dir):
        with warnings.catch_warnings():  # turning off ResourceWarnings (unclosed file) from ete3
            warnings.simplefilter('ignore')
            tax_path = NCBITaxonomyDb._define_tax_paths(data_dir)
            if not os.path.exists(tax_path):
                logging.error('NCBI files not found in specified directory. Use metaquantome db to download files')
            return NCBITaxa(tax_path)

    def is_in_db(self, taxid):
        # todo: doc
        rank_dict = self.ncbi.get_rank([taxid])
        if taxid in rank_dict.keys():
            return True
        else:
            return False

    # @staticmethod
    # def _ncbi_database_handler(data_dir):
    #     # todo: doc
    #     with warnings.catch_warnings():  # turning off ResourceWarnings (unclosed file) from ete3
    #         warnings.simplefilter('ignore')
    #         tax_path = os.path.join(data_dir, 'taxa.sqlite')
    #         if os.path.exists(tax_path):
    #             logging.info('Using taxonomy database in ' + tax_path)
    #             ncbi = NCBITaxa(tax_path)
    #         else:
    #             open(tax_path, 'a').close()  # will say that database is not up to date
    #             ncbi = NCBITaxa(tax_path)  # this prints a lot of logging messages, so no message here
    #             # remove taxdump, in working directory
    #             taxdump = os.path.join(os.getcwd(), 'taxdump.tar.gz')
    #             if os.path.exists(taxdump):
    #                 os.remove(taxdump)
    #     return ncbi

    def map_id_to_desired_ranks(self, ranks2get, taxid):
        # todo: doc
        rank_of_query = self.get_rank(taxid)
        num_rank_of_query = NUMERIC_RANK[rank_of_query]
        lineage = self.ncbi.get_lineage(taxid)
        ranks = self.ncbi.get_rank(lineage)
        # flip the dict so we have rank: id, rather than id: rank
        # assumes we have unique values for all relevant ranks (i.e., order)
        # which should be true
        invert_dict = {v: k for k, v in ranks.items()}
        # get ranks and taxids
        # if rank is not present in lineage, say is unidentified
        lineage_ranks = invert_dict.keys()
        rt_dict = dict()
        for rank in ranks2get:
            num_rank_of_rank2get = NUMERIC_RANK[rank]
            # must be above taxid in hierarchy
            if num_rank_of_rank2get <= num_rank_of_query:
                if rank in lineage_ranks:
                    rt_dict[rank] = invert_dict[rank]
                else:
                    rt_dict[rank] = UNIDENTIFIED
        return rt_dict

    def expand_sample_taxonomy(self, sample_set):
        """
        Expand sample set to include ancestors
        :param sample_set: set of all taxids explicitly provided
        :return: set of all taxids, with all ancestors
        """
        # base set
        expanded_sample_set = sample_set.copy()
        # for each term, get all (relevant) ranks above it.
        for i in sample_set:
            expanded_sample_set.update(set(self.map_id_to_desired_ranks(BASIC_TAXONOMY_TREE, i).values()))
        return expanded_sample_set

    def filter_to_desired_ranks(self, taxids, desired_ranks=BASIC_TAXONOMY_TREE):
        # todo: doc
        ranks = self.ncbi.get_rank(taxids)
        relevant_ranks = {k for k, v in ranks.items() if v in desired_ranks}
        return relevant_ranks

    def get_rank(self, taxid):
        # todo: doc
        # get the rank of the query, as a number
        query_rank = self.ncbi.get_rank([taxid])[taxid]
        return query_rank

    def get_children(self, taxid):
        # todo: doc
        # get the rank of the query, as a number
        query_rank = self.get_rank(taxid)
        # get number of query rank
        num_query_rank = NUMERIC_RANK[query_rank]
        # which major rank is the first with lower rank?
        try:
            child_rank_char = BASIC_TAXONOMY_TREE[[v > num_query_rank for v in BASIC_NUMERIC_RANK].index(True)]
        except ValueError:  # if rank is species (none lower in BASIC_NUMERIC_RANK)
            return set()
        # get all descendants of query
        desc = self.ncbi.get_descendant_taxa(taxid, intermediate_nodes=True)
        # get ranks of all descendants
        desc_ranks = self.ncbi.get_rank(desc)
        # filter all descendants by rank
        immed_children = {k for k, v in desc_ranks.items() if v == child_rank_char}
        return immed_children

    def get_descendants(self, taxid):
        # todo: doc
        all_descendants = self.ncbi.get_descendant_taxa(taxid, intermediate_nodes=True)
        # only return those in BASIC_TAXONOMY_TREE
        relevant_descendants = self.filter_to_desired_ranks(all_descendants, desired_ranks=BASIC_TAXONOMY_TREE)
        return relevant_descendants

    def get_parents(self, taxid):
        # todo: doc
        rank_of_query = self.get_rank(taxid)  # needs to be list
        num_query_rank = NUMERIC_RANK[rank_of_query]
        # which major rank is the lowest with lower rank? (i.e., more general)
        # we want the *largest* numerical rank that is above the query rank
        try:
            lower_list = [v < num_query_rank for v in BASIC_NUMERIC_RANK]
            rev_lower_list = lower_list[::-1]  # reverse, so the first occurrence of True is now the lowest rank
            # get the index of the first occurrence of True, then subtract it from BASIC_TAXONOMY_TREE
            parent_rank = len(BASIC_TAXONOMY_TREE) - rev_lower_list.index(True) - 1
            parent_rank_char = BASIC_TAXONOMY_TREE[parent_rank]
        except ValueError:  # if rank is phylum (none higher in BASIC_NUMERIC_RANK)
            return set()
        # get full lineage of query
        lineage = self.ncbi.get_lineage(taxid)
        # get only parents
        parents = self.filter_to_desired_ranks(lineage, [parent_rank_char])
        return parents

    def get_ancestors(self, taxid):
        # todo: doc
        rank_of_query = self.get_rank(taxid)
        num_query_rank = NUMERIC_RANK[rank_of_query]
        # which ranks are more general?
        try:
            ancestor_rank_char = [BASIC_TAXONOMY_TREE[i] for
                                  i, v in enumerate(BASIC_NUMERIC_RANK) if v < num_query_rank]
        except ValueError:  # if rank is phylum (none higher in BASIC_NUMERIC_RANK)
            return set()
        # get full lineage of query
        lineage = self.ncbi.get_lineage(taxid)
        # get only ancestors
        ancestors = self.filter_to_desired_ranks(lineage, ancestor_rank_char)
        return ancestors

    def convert_taxid_to_name(self, taxids):
        # todo: doc
        """
        :param taxids: a list or Series of taxids
        :return: list of ncbi names
        """
        translator = self.ncbi.get_taxid_translator(taxids)
        # id must be an integer, so cast from string
        names = [translator[int(tid)] for tid in taxids]
        return names

    def convert_name_to_taxid(self, names):
        # todo: doc
        """
        Converts the name of a given taxon to its numerical NCBI taxonomy id.
        If a name maps to multiple IDs, the first in the list is taken - this is somewhat unpredictable,
        so using taxonomic id's is strongly preferred.
        If the name is not found in the database, than it is replaced with NaN, which will
        be filtered out later.
        :param names: a list or taxon names (as strings)
        :return: list of ncbi taxonomy ids
        """
        # arg checking
        if not isinstance(names, list):
            TypeError('names must be list')

        translator = self.ncbi.get_name_translator(names)
        ids = [0] * len(names)
        for i in range(len(names)):
            if names[i] in translator.keys():
                ids[i] = translator[names[i]][0]  # always takes first id
            else:
                ids[i] = np.nan  # make missing, will be filtered out later
        return ids
