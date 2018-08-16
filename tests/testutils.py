import os
from metaquant.definitions import DATA_DIR
import pandas

pandas.set_option('display.max_columns', 20)

def testfile(name):
    return os.path.join(DATA_DIR, 'test', name)

