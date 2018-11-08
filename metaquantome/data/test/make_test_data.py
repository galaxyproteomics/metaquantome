import pandas as pd
import os
from metaquantome.util.utils import DATA_DIR
import numpy as np


def write_testfile(df, name):
    df.to_csv(os.path.join(DATA_DIR, 'test', name), sep='\t', index_label='peptide')


# simple: single intensity

func = pd.DataFrame({'go': ['GO:0008152', 'GO:0022610']}, index=['A', 'B'])
write_testfile(func, 'simple_func.tab')

ec = pd.DataFrame({'ec': ['3.4.11.-', '3.4.21.70']}, index=['A', 'B'])
write_testfile(ec, 'simple_ec.tab')

unk_ec = pd.DataFrame({'ec': ['1.50.10000.-', '3.4.21.70']}, index=['A', 'B'])
write_testfile(unk_ec, 'unk_ec.tab')

int = pd.DataFrame({'int': ['100', '200']}, index=['A', 'B'])
write_testfile(int, 'simple_int.tab')

# taxa are h. pylori and c. difficile (both species), from different phyla
tax = pd.DataFrame({'lca': [210, 1496]}, index=['A', 'B'])
write_testfile(tax, 'simple_tax.tab')

# multiple intensities
peptides = ['A', 'B', 'C']

mult_func = pd.DataFrame({'go': ['GO:0008152', 'GO:0022610', 'GO:0000003,GO:0032505'],
                          'cog': ['C', 'N', 'D'],
                          'ec': ['3.4.11.-', '3.4.21.70', '1.2.-.-']},
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
t_int = pd.DataFrame({'int1': [12, 20, 1000],
                      'int2': [20, 30, 1200],
                      'int3': [15, 20, 900],
                      'int4': [12, 3500, 12],
                      'int5': [21, 2000, 13],
                      'int6': [10, 3000, 10]}, index=peptides)
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

# clean up unipept_sample7_functional.csv
unipept_direct='unipept_sample7_functional.csv'
uni_raw = pd.read_csv(unipept_direct, index_col='peptide')
uni = uni_raw[['EC']].copy()
# converts sep within column to ',', removes parentheses
uni.replace({'EC': {' \(\d{1,3}\%\)': '', ';':','}}, regex=True, inplace=True)

# fake intensities
np.random.seed(101)
uni['int'] = np.random.lognormal(10, 1, uni.size)
# set anything with 1 as top level category to 0, for testing
uni.loc[uni.EC.str.contains(r'^1|,1',na=False), 'int'] = 0

uni.to_csv(os.path.join(DATA_DIR, 'test', 'unipept_sample7_functional_clean.tab'), sep='\t', columns=['EC'],
           index_label="peptide")
uni.to_csv(os.path.join(DATA_DIR, 'test', 'unipept_sample7_int_clean.tab'), sep='\t', columns=['int'],
           index_label="peptide")

# make lca table for testing filtering
uni_raw.to_csv(os.path.join(DATA_DIR, 'test', 'unipept_sample7_taxonomy.tab'), sep='\t', columns=['lca'],
           index_label="peptide")

# make nopep table

nopep = pd.DataFrame({'go': ['GO:0008152', 'GO:0022610', 'GO:0000003,GO:0032505'],
                      'cog': ['C', 'N', None,],
                      'ec': ['3.4.11.-', '3.4.21.70', '1.2.-.-'],
                      'int1': [10, 40, 50],
                      'int2': [20, 30, 50],
                      'int3': [70, 30, 0],
                      'lca': ['210', '1496', '1870884']})
nopep.to_csv(os.path.join(DATA_DIR, 'test', 'nopep.tab'), sep='\t',
             index=False)


# for filtering

# h pylori, c difficile, and clostridiaceae family
# t-testing: 6 samples, intensities very different
filt_int = pd.DataFrame({'int1': [12, 20, 1000],
                         'int2': [20, 0, 0],
                         'int3': [15, 20, 900],
                         'int4': [12, 0, 12],
                         'int5': [21, 2000, 13],
                         'int6': [10, 3000, 10]}, index=peptides)
write_testfile(filt_int, 'filt_int.tab')
