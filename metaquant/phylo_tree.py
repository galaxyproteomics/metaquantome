from ete3 import NCBITaxa
import os
import numpy as np


def load_ncbi(update=False):
    ncbi = NCBITaxa()
    if update:
        # database is saved to ~/.etetoolkit / taxa.sqlite
        ncbi.update_taxonomy_database()

        # remove taxdump, in working directory
        taxdump = os.path.join(os.getcwd(), 'taxdump.tar.gz')
        if os.path.exists(taxdump):
            os.remove(taxdump)
            # if not found, do nothing
    return ncbi


def get_desired_ranks_from_lineage(ranks2get, taxid, ncbi):

    taxid = handle_nan_taxid(taxid)

    lineage = ncbi.get_lineage(taxid)
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

