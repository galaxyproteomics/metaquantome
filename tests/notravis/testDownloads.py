import unittest
import os
import shutil

from metaquantome.databases.EnzymeDb import EnzymeDb
from metaquantome.databases.GeneOntologyDb import GeneOntologyDb
import metaquantome.databases.NCBITaxonomyDb as td
from metaquantome.util.utils import DATA_DIR


class TestDownloads(unittest.TestCase):
    def testEnzymeDatabaseHandler(self):
        # make tmp dir for testing download
        tmp_dir = os.path.join(DATA_DIR, 'tmp_test_data_dwnld')
        os.mkdir(tmp_dir)
        try:
            EnzymeDb.download_enzyme_db(tmp_dir, False)
            enzyme_db = EnzymeDb(tmp_dir)
            expected_contents = [os.path.join(tmp_dir, file)
                                 for file in ['enzclass.txt', 'enzyme.dat', 'ec_id.json', 'enzclass.json']]
            for content in expected_contents:
                self.assertTrue(os.path.exists(content))
            # make sure parsed correctly
            # this is from enzyme.dat
            self.assertEqual(enzyme_db.ecdb['1.2.3.4']['descript'], 'Oxalate oxidase.')

            # from enzclass.txt
            self.assertEqual(enzyme_db.ecdb['6.1.-.-']['descript'], 'Forming carbon-oxygen bonds.')
        finally:
            shutil.rmtree(tmp_dir)

    def testDownloadGO(self):
        tmp_dir = os.path.join(DATA_DIR, 'tmp_go_data_dwnld')
        os.mkdir(tmp_dir)
        try:
            GeneOntologyDb.download_go(tmp_dir, overwrite=False)
            go = GeneOntologyDb(tmp_dir, slim_down=True)
            expected_contents = [os.path.join(tmp_dir, file)
                                 for file in ['go-basic.obo', 'goslim_metagenomics.obo']]
            for content in expected_contents:
                self.assertTrue(os.path.exists(content))
            # make sure parsed correctly
            self.assertEqual('biological_process', go.gofull['GO:0008150'].name)
        finally:
            shutil.rmtree(tmp_dir)

    def testDownloadTaxonomy(self):
        # make tmp dir for testing download
        tmp_dir = os.path.join(DATA_DIR, 'tmp_test_tax_dwnld')
        os.mkdir(tmp_dir)
        try:
            td.NCBITaxonomyDb.download_ncbi(tmp_dir)
            ncbi = td.NCBITaxonomyDb(tmp_dir)
            lineage = ncbi.get_ancestors(1919)
            self.assertTrue(1760 in lineage)
        finally:
            shutil.rmtree(tmp_dir)
