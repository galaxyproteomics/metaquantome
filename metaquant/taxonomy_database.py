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


class NCBITaxonomyDb:
    # todo: implement get all ancestors

    def __init__(self, data_dir):
        self.ncbi = self._ncbi_database_handler(data_dir)

    def handle_nan_taxid(self, taxid):
        # replace nan's with 'unidentified'
        # only nan if not character
        if isinstance(taxid, str):
            taxid = int(taxid)

        if np.isnan(taxid):
            taxid = 32644

        return taxid

    def _ncbi_database_handler(self, data_dir):
        with warnings.catch_warnings():  # turning off ResourceWarnings (unclosed file) from ete3
            warnings.simplefilter('ignore')
            ncbi_db_path = os.path.join(data_dir, 'ncbi')
            if not os.path.exists(ncbi_db_path):
                os.mkdir(ncbi_db_path)
            tax_path = os.path.join(ncbi_db_path, 'taxa.sqlite')
            if os.path.exists(tax_path):
                logging.info('Using taxonomy database in ' + tax_path)
                ncbi = NCBITaxa(tax_path)
            else:
                open(tax_path, 'a').close()  # will say that database is not up to date
                ncbi = NCBITaxa(tax_path)  # this prints a lot of logging messages, so no message here
                # remove taxdump, in working directory
                taxdump = os.path.join(os.getcwd(), 'taxdump.tar.gz')
                if os.path.exists(taxdump):
                    os.remove(taxdump)
        return ncbi

    def map_id_to_desired_ranks(self, ranks2get, taxid):

        clean_taxid = self.handle_nan_taxid(taxid)
        rank_of_query = self.ncbi.get_rank([clean_taxid])[clean_taxid]  # needs to be list
        num_rank_of_query = NUMERIC_RANK[rank_of_query]

        lineage = self.ncbi.get_lineage(clean_taxid)
        ranks = self.ncbi.get_rank(lineage)

        # flip the dict so we have rank: id, rather than id: rank
        # assumes we have unique values for all relevant ranks (i.e., order)
        # which should be true
        invert_dict = {v: k for k, v in ranks.items()}

        # get ranks and taxids
        # if rank is not present in lineage, say is unclassified
        lineage_ranks = invert_dict.keys()
        rt_dict = dict()
        for rank in ranks2get:
            num_rank_of_rank2get = NUMERIC_RANK[rank]
            # must be above taxid in hierarchy
            if num_rank_of_rank2get <= num_rank_of_query:
                if rank in lineage_ranks:
                    rt_dict[rank] = invert_dict[rank]
                else:
                    rt_dict[rank] = 12908  # ncbi taxid for unclassified
        return rt_dict

    def expand_sample_taxonomy(self, sample_set):
        """
        Expand sample set to include ancestors
        :param sample_set: set of all taxids explicitly provided
        :param ncbi: NCBITaxa() object
        :return: set of all taxids, with all ancestors
        """

        # base set
        expanded_sample_set = sample_set.copy()

        # for each term, get all (relevant) ranks above it.
        for i in sample_set:
            expanded_sample_set.update(set(self.map_id_to_desired_ranks(BASIC_TAXONOMY_TREE, i).values()))

        return expanded_sample_set

    def filter_to_desired_ranks(self, taxids, desired_ranks=BASIC_TAXONOMY_TREE):
        ranks = self.ncbi.get_rank(taxids)
        relevant_ranks = {k for k, v in ranks.items() if v in desired_ranks}
        return relevant_ranks

    def get_children(self, id):
        # get the rank of the query, as a number
        query_rank = self.ncbi.get_rank([id])[id]

        # get number of query rank
        num_query_rank = NUMERIC_RANK[query_rank]

        # which major rank is the first with lower rank?
        try:
            child_rank_char = BASIC_TAXONOMY_TREE[[v > num_query_rank for v in BASIC_NUMERIC_RANK].index(True)]
        except ValueError:  # if rank is species (none lower in BASIC_NUMERIC_RANK)
            return set()

        # get all descendants of query
        desc = self.ncbi.get_descendant_taxa(id, intermediate_nodes=True)

        # get ranks of all descendants
        desc_ranks = self.ncbi.get_rank(desc)

        # filter all descendants by rank
        immed_children = {k for k,v in desc_ranks.items() if v == child_rank_char}
        return immed_children

    def get_descendants(self, id):
        all_descendants = self.ncbi.get_descendant_taxa(id, intermediate_nodes=True)

        # only return those in BASIC_TAXONOMY_TREE
        relevant_descendants = self.filter_to_desired_ranks(all_descendants, desired_ranks=BASIC_TAXONOMY_TREE)
        return relevant_descendants

    def get_parents(self, id):
        clean_taxid = self.handle_nan_taxid(id)
        rank_of_query = self.ncbi.get_rank([clean_taxid])[clean_taxid]  # needs to be list
        num_query_rank = NUMERIC_RANK[rank_of_query]

        # which major rank is the lowest with lower rank? (i.e., more general)
        # we want the *largest* numerical rank that is above the query rank
        try:
            lower_list = [v < num_query_rank for v in BASIC_NUMERIC_RANK]
            rev_lower_list = lower_list[::-1] # reverse, so the first occurence of True is now the lowest rank
            # get the index of the first occurence of True, then subtract it from BASIC_TAXONOMY_TREE
            parent_rank = len(BASIC_TAXONOMY_TREE) - rev_lower_list.index(True) - 1
            parent_rank_char = BASIC_TAXONOMY_TREE[parent_rank]
        except ValueError:  # if rank is phylum (none higher in BASIC_NUMERIC_RANK)
            return set()

        # get full lineage of query
        lineage = self.ncbi.get_lineage(id)

        # get only parents
        parents = self.filter_to_desired_ranks(lineage, [parent_rank_char])
        return parents

    def get_ancestors(self, id):
        clean_taxid = self.handle_nan_taxid(id)
        rank_of_query = self.ncbi.get_rank([clean_taxid])[clean_taxid]  # needs to be list
        num_query_rank = NUMERIC_RANK[rank_of_query]

        # which ranks are more general?
        try:
            ancestor_rank_char = [BASIC_TAXONOMY_TREE[i] for
                                  i,v in enumerate(BASIC_NUMERIC_RANK) if v < num_query_rank]
        except ValueError:  # if rank is phylum (none higher in BASIC_NUMERIC_RANK)
            return set()

        # get full lineage of query
        lineage = self.ncbi.get_lineage(id)

        # get only ancestors
        ancestors = self.filter_to_desired_ranks(lineage, ancestor_rank_char)
        return ancestors

    def convert_taxid_to_name(self, taxids):
        """
        :param taxids: a list or Series of taxids
        :param ncbi: a NCBITaxa object, from ete3
        :return: list of ncbi names
        """
        # fill nans with unassigned taxid
        safe_taxids = [self.handle_nan_taxid(id) for id in taxids]

        translator = self.ncbi.get_taxid_translator(safe_taxids)

        # id must be an integer, so cast from string
        names = [translator[int(id)] for id in safe_taxids]
        return names

    def convert_name_to_taxid(self, names):
        """

        :param names: a list or series of taxon names (as strings)
        :param ncbi: a NCBITaxa object from ete3 package
        :return: list of ncbi taxonomy ids
        """
        translator = self.ncbi.get_name_translator(names)
        ids = [0] * len(names)

        for i in range(len(names)):
            if names[i] in translator.keys():
                ids[i] = translator[names[i]][0]  # always takes first id
            else:
                ids[i] = 32644  # ncbi taxid for unidentified

        return ids


# the following two are for all 4 database types
# todo: move to sample set file
def number_of_children(id, db, dbtype, expanded_sample_set):
    """
    get number of children within sample set
    :param id:
    :param db:
    :param dbtype:
    :param expanded_sample_set: set of all taxa in sample, with all ancestors
    :return: the number of children for the id node
    """
    if dbtype == 'ncbi':
        int_taxid = int(id)
        if id in {1, 2, 2759, 2157}:  # root, bacteria, eukaryota, or archaea
            return np.inf
        else:
            children = db.get_children(id)
            nchildren = len(children.intersection(expanded_sample_set))
            return nchildren


def prune_taxonomy_tree(expanded_sample_set, min_children, ncbi):
    # nchildren
    good_ids = expanded_sample_set.copy()
    for i in expanded_sample_set:
        nchild = number_of_children(i, ncbi, 'ncbi', expanded_sample_set)
        if 0 < nchild < min_children:
            good_ids.remove(i)
    return good_ids
