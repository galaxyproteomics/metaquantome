import os

from metaquantome.util.utils import DATA_DIR

MISSING_VALUES = ["", "0", "NA", "NaN", "0.0"]
ONTOLOGIES = ['cog', 'go', 'ec']
P_COLNAME = 'p'
P_CORR_COLNAME = 'corrected_p'
TAX_TEST_DIR = os.path.join(DATA_DIR, 'test', 'ncbi_test')
GO_TEST_DIR = os.path.join(DATA_DIR, 'test', 'go_cache')  # downloaded 11/5/18
EC_TEST_DIR = os.path.join(DATA_DIR, 'test', 'ec_cache')  # downloaded 8/28/18
