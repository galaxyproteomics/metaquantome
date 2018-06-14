from ete3 import NCBITaxa
import os


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

    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)

    # assumes we have unique values for all relevant ranks (i.e., order)
    # which should be true
    invert_dict = {v : k for k, v in ranks.items()}

    # get ranks and taxids
    # if rank is not present in lineage, say is 'unknown'
    lineage_ranks = invert_dict.keys()
    rt_dict = dict()
    for rank in ranks2get:
        if rank in lineage_ranks:
            rt_dict[rank] = invert_dict[rank]
        else:
            rt_dict[rank] = 32644  # ncbi taxid for unassigned

    return rt_dict


def convert_taxid_to_name(taxid, ncbi):
    translator = ncbi.get_taxid_translator(taxid)

    # id must be an integer, so cast
    names = [translator[int(id)] for id in taxid]
    return names


def convert_name_to_taxid(names, ncbi):
    translator = ncbi.get_name_translator(names)
    ids = [0]*len(names)

    for i in range(len(names)):
        if names[i] in translator.keys():
            ids[i] = translator[names[i]][0]  # always takes first id
        else:
            ids[i] = 32644  # ncbi taxid for unassigned

    return ids

