# todo: add check args
"""
check that hard tabs are present in all tabular files
make sure database files are present in a more elegant way
make sure text strings match potential options
others?
"""
import os


def function_check(func_file, func_colname):
    # todo: doc
    if not func_file:
        raise IOError('Function tabular file not provided (--func_file)')
    if not os.path.exists(func_file):
        raise FileNotFoundError('func_file does not exist. Please check filename')
    if not func_colname:
        raise ValueError('func_colname=None. Please provide a function column name (--func_colname)')


def tax_check(tax_file, tax_colname):
    # todo: doc
    if not tax_file:
        raise IOError('Taxonomy tabular file not provided (--tax_file)')
    if not os.path.exists(tax_file):
        raise FileNotFoundError('func_file does not exist. Please check filename')
    if not tax_colname:
        raise ValueError('tax_colname=None. Please provide a taxonomy column name (--tax_colname)')