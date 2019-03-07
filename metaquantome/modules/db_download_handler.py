from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
from metaquantome.databases.EnzymeDb import EnzymeDb
from metaquantome.databases.NCBITaxonomyDb import NCBITaxonomyDb


def db_download_handler(dbs, data_dir, overwrite):
    """
    called by CLI to download databases

    :param dbs: list of databases ("go", "ec", and/or "ncbi")
    :param data_dir: data directory
    :param overwrite: whether to overwrite existing files or not
    :return: None
    """
    if "go" in dbs:
        GeneOntologyDb.download_go(data_dir, overwrite)
    if "ec" in dbs:
        EnzymeDb.download_enzyme_db(data_dir, overwrite)
    if "ncbi" in dbs:
        NCBITaxonomyDb.download_ncbi(data_dir)
