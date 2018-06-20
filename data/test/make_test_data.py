import pandas as pd
import os
from definitions import DATA_DIR


def write_testfile(df, name):
    df.to_csv(os.path.join(DATA_DIR, 'test', name), sep='\t', index_label='peptide')


# simple: single intensity

func = pd.DataFrame({'go': ['GO:0008152', 'GO:0022610']}, index=['A', 'B'])
write_testfile(func, 'simple_func.tab')

int = pd.DataFrame({'int': ['100', '200']}, index=['A', 'B'])
write_testfile(int, 'simple_int.tab')

# taxa are h. pylori and c. difficile
tax = pd.DataFrame({'lca': ['210', '1496']}, index=['A', 'B'])
write_testfile(tax, 'simple_tax.tab')

# multiple intensities
peptides = ['A', 'B', 'C']

mult_func = pd.DataFrame({'go': ['GO:0008152', 'GO:0022610', 'GO:0000003,GO:0032505'],
                          'cog': ['C', 'N', 'D']},
                         index=peptides)
write_testfile(mult_func, 'multiple_func.tab')

mult_int = pd.DataFrame({'int1': [10, 40, 50],
                         'int2': [20, 30, 50],
                         'int3': [70, 30, 0]},
                        index=peptides)
write_testfile(mult_int, 'multiple_int.tab')

# h pylori, c difficile, and clostridiaceae family
mult_tax = pd.DataFrame({'lca': ['210', '1496', '1870884']}, index=peptides)
write_testfile(mult_tax, 'multiple_tax.tab')

# t-testing: 6 samples, intensities very different
t_int = pd.DataFrame({'int1': [10, 20, 1000],
                      'int2': [20, 30, 1200],
                      'int3': [15, 20, 900],
                      'int4': [30, 3500, 12],
                      'int5': [21, 2000, 13],
                      'int6': [30, 3000, 10]}, index=peptides)
write_testfile(t_int, 'int_ttest.tab')


# funtax interaction

mult_ft_func = pd.DataFrame({'cog': ['C', 'N', 'C']}, index=peptides)
write_testfile(mult_ft_func, 'mult_ft_func.tab')

mult_ft_tax = pd.DataFrame({'lca': ['210', '1496', '210']}, index=peptides)
write_testfile(mult_ft_tax, 'mult_ft_tax.tab')

# sample information file
samp_info = pd.DataFrame({'group': ['A', 'B'],
                          'colnames': ['A1,A2', 'B1, B2'],
                          })

samp_info.to_csv(os.path.join(DATA_DIR, 'test', 'samp_info.tab'),
                 sep='\t',
                 columns=['group', 'colnames'], index=False)
