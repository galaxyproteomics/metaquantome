import os
from metaquant.definitions import DATA_DIR


def testfile(name):
    return os.path.join(DATA_DIR, 'test', name)

