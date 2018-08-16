from ete3 import NCBITaxa
import os
import numpy as np
import logging
import warnings


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


def ncbi_database_handler(data_dir):
    with warnings.catch_warnings(): # turning off ResourceWarnings (unclosed file) from ete3
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


def map_id_to_desired_ranks(ranks2get, taxid, ncbi):

    clean_taxid = handle_nan_taxid(taxid)
    rank_of_query = ncbi.get_rank([clean_taxid])[clean_taxid] # needs to be list
    num_rank_of_query = NUMERIC_RANK[rank_of_query]

    lineage = ncbi.get_lineage(clean_taxid)
    ranks = ncbi.get_rank(lineage)

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


def handle_nan_taxid(taxid):
    # replace nan's with 'unidentified'
    # only nan if not character
    if isinstance(taxid, str):
        taxid = int(taxid)

    if np.isnan(taxid):
        taxid = 32644

    return taxid


def convert_taxid_to_name(taxids, ncbi):
    """
    :param taxids: a list or Series of taxids
    :param ncbi: a NCBITaxa object, from ete3
    :return: list of ncbi names
    """
    # fill nans with unassigned taxid
    safe_taxids = [handle_nan_taxid(id) for id in taxids]

    translator = ncbi.get_taxid_translator(safe_taxids)

    # id must be an integer, so cast from string
    names = [translator[int(id)] for id in safe_taxids]
    return names


def convert_name_to_taxid(names, ncbi):
    """

    :param names: a list or series of taxon names (as strings)
    :param ncbi: a NCBITaxa object from ete3 package
    :return: list of ncbi taxonomy ids
    """
    translator = ncbi.get_name_translator(names)
    ids = [0]*len(names)

    for i in range(len(names)):
        if names[i] in translator.keys():
            ids[i] = translator[names[i]][0]  # always takes first id
        else:
            ids[i] = 12908  # ncbi taxid for unclassified

    return ids


def expand_sample_taxonomy(sample_set, ncbi):
    '''
    Expand sample set to include ancestors
    :param sample_set: set of all taxids explicitly provided
    :param ncbi: NCBITaxa() object
    :return: set of all taxids, with all ancestors
    '''
    expanded_sample_set = sample_set

    # for each term, get all ranks above it.

