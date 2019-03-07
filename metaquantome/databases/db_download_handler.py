from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
from metaquantome.databases.EnzymeDb import EnzymeDb
from metaquantome.databases.NCBITaxonomyDb import NCBITaxonomyDb


def db_download_handler(dbs, dir, overwrite):
    if "go" in dbs:
        GeneOntologyDb.download_go(dir, overwrite)
    if "ec" in dbs:
        EnzymeDb.download_enzyme_db(dir, overwrite)
    if "ncbi" in dbs:
        NCBITaxonomyDb.download_ncbi(dir)
