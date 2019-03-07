import os

from metaquantome.util.utils import DATA_DIR

MISSING_VALUES = ["", "0", "NA", "NaN", "0.0"]
ONTOLOGIES = ['cog', 'go', 'ec']
P_COLNAME = 'p'
P_CORR_COLNAME = 'corrected_p'

# go database downloaded 11/5/18
# ec database downloaded 8/28/18
TEST_DIR = os.path.join(DATA_DIR, 'test')

