import os
from metaquant.utils import DATA_DIR
import pandas

pandas.set_option('display.max_columns', 20)


def testfile(name):
    return os.path.join(DATA_DIR, 'test', name)


TTEST_SINFO = '{"s1": ["int1", "int2", "int3"], "s2": ["int4", "int5", "int6"]}'
