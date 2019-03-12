import os


def function_check(func_file, func_colname):
    """
    check that the function file is provided and it exists.
    check that the function colname is provided

    :param func_file: path to function file
    :param func_colname: string - name of function column in func_file
    :return: None
    """
    if not func_file:
        raise IOError('Function tabular file not provided (--func_file)')
    if not os.path.exists(func_file):
        raise FileNotFoundError('func_file does not exist. Please check filename')
    if not func_colname:
        raise ValueError('func_colname=None. Please provide a function column name (--func_colname)')


def tax_check(tax_file, tax_colname):
    """
    check that the taxonomy file is provided and it exists.
    check that the taxonomy colname is provided

    :param tax_file: path to taxonomy file
    :param tax_colname: string - name of taxonomy column in tax_file
    :return: None
    """
    if not tax_file:
        raise IOError('Taxonomy tabular file not provided (--tax_file)')
    if not os.path.exists(tax_file):
        raise FileNotFoundError('func_file does not exist. Please check filename')
    if not tax_colname:
        raise ValueError('tax_colname=None. Please provide a taxonomy column name (--tax_colname)')
