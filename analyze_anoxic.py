import metaquant

# for anoxic lake, don't do tests because we don't have different conditions
# more exploratory analysis

infile = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/anoxic_lake/anoxic_all_combined.tabular'
#
# fun1outfile = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/anoxic_lake/anoxic_function1_out.tabular'
# fun2outfile = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/anoxic_lake/anoxic_function2_out.tabular'
#
# f = metaquant.metaquant('fn',
#                         file=infile,
#                         sample1_colnames=['sz9_1', 'sz9_2', 'sz9_3'],
#                         sample2_colnames=['sz10_1', 'sz10_2', 'sz10_3'],
#                         func_colname='go',
#                         pep_colname='peptide',
#                         ontology="GO",
#                         test=False,
#                         outfile=fun1outfile)
#
# f2 = metaquant.metaquant('fn',
#                         file=infile,
#                         sample1_colnames=['sz9_1', 'sz9_2', 'sz9_3'],
#                         sample2_colnames=['sz11_1', 'sz11_2', 'sz11_3'],
#                         func_colname='go',
#                         pep_colname='peptide',
#                         ontology="GO",
#                         test=False,
#                         outfile=fun2outfile)

# funslim1_outfile = '/home/caleb/Griffin_lab_work/' +\
#                   'functional_taxonomic_quant/anoxic_lake/anoxic_lake_function_out1_goslim.tabular'
# funslim2_outfile = '/home/caleb/Griffin_lab_work/' +\
#                   'functional_taxonomic_quant/anoxic_lake/anoxic_lake_function_out2_goslim.tabular'
#
# fslim1 = metaquant.metaquant('fn',
#                             file=infile,
#                             sample1_colnames=['sz9_1', 'sz9_2', 'sz9_3'],
#                             sample2_colnames=['sz10_1', 'sz10_2', 'sz10_3'],
#                             func_colname='go',
#                             pep_colname='peptide',
#                             ontology="GO",
#                             slim_down=True,
#                             test=False,
#                             outfile=funslim1_outfile)
#
# fslim2 = metaquant.metaquant('fn',
#                         file=infile,
#                         sample1_colnames=['sz9_1', 'sz9_2', 'sz9_3'],
#                         sample2_colnames=['sz11_1', 'sz11_2', 'sz11_3'],
#                         func_colname='go',
#                         pep_colname='peptide',
#                         ontology="GO",
#                         slim_down=True,
#                         test=False,
#                         outfile=funslim2_outfile)
# #
# # cog
# fcog_outfile = '/home/caleb/Griffin_lab_work/' +\
#                   'functional_taxonomic_quant/anoxic_lake/anoxic_lake_function_out_cog.tabular'
# fcog = metaquant.metaquant('fn',
#                             file=infile,
#                             sample1_colnames=['sz9_1', 'sz9_2', 'sz9_3'],
#                             sample2_colnames=['sz10_1', 'sz10_2', 'sz10_3'],
#                             func_colname='cog',
#                             pep_colname='peptide',
#                             ontology="cog",
#                             test=False,
#                             outfile=fcog_outfile)
#
taxout1 = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/anoxic_lake/anoxic_lake_taxonomy_maxquant_out1.tabular'
taxout2 = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/anoxic_lake/anoxic_lake_taxonomy_maxquant_out2.tabular'

t1 = metaquant.metaquant('tax',
                        file=infile,
                        sample1_colnames=['sz9_1', 'sz9_2', 'sz9_3'],
                        sample2_colnames=['sz10_1', 'sz10_2', 'sz10_3'],
                        outfile=taxout1)

t2 = metaquant.metaquant('tax',
                        file=infile,
                        sample1_colnames=['sz9_1', 'sz9_2', 'sz9_3'],
                        sample2_colnames=['sz11_1', 'sz11_2', 'sz11_3'],
                        outfile=taxout2)

#
# # fun tax interaction
# ftout = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/anoxic_lake/anoxic_lake_funtax_maxquant_out.tabular'
# ft = metaquant.metaquant('taxfn',
#                          file=infile,
#                          sample1_colnames=['sz9_1', 'sz9_2', 'sz9_3'],
#                          sample2_colnames=['sz10_1', 'sz10_2', 'sz10_3'],
#                          lca_colname="lca",
#                          cog_colname="cog",
#                          outfile=ftout)