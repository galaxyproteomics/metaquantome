import metaquant


infile = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/rudney3pairs/function_taxonomy_intensity_maxquant_unipept.tabular'

# funoutfile = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/rudney3pairs/rudney_3pairs_function_out.tabular'
#
# f = metaquant.metaquant('fn',
#                         file=infile,
#                         sample1_colnames=['Intensity.737NS', 'Intensity.852NS', 'Intensity.867NS'],
#                         sample2_colnames=['Intensity.737WS', 'Intensity.852WS', 'Intensity.867WS'],
#                         func_colname='go',
#                         pep_colname='peptide',
#                         ontology="GO",
#                         test=True,
#                         paired=True,
#                         outfile=funoutfile,
#                         threshold=3)
#
# funslim_outfile = '/home/caleb/Griffin_lab_work/' +\
#                   'functional_taxonomic_quant/rudney3pairs/rudney_3pairs_function_out_goslim.tabular'
#
# fslim = metaquant.metaquant('fn',
#                             file=infile,
#                             sample1_colnames=['Intensity.737NS', 'Intensity.852NS', 'Intensity.867NS'],
#                             sample2_colnames=['Intensity.737WS', 'Intensity.852WS', 'Intensity.867WS'],
#                             func_colname='go',
#                             pep_colname='peptide',
#                             ontology="GO",
#                             slim_down=True,
#                             test=True,
#                             paired=True,
#                             outfile=funslim_outfile,
#                             threshold=3)

# cog
fcog_outfile = '/home/caleb/Griffin_lab_work/' +\
                  'functional_taxonomic_quant/rudney3pairs/rudney_3pairs_function_out_cog.tabular'
fcog = metaquant.metaquant('fn',
                            file=infile,
                            sample1_colnames=['Intensity.737NS', 'Intensity.852NS', 'Intensity.867NS'],
                            sample2_colnames=['Intensity.737WS', 'Intensity.852WS', 'Intensity.867WS'],
                            func_colname='cog',
                            pep_colname='peptide',
                            ontology="cog",
                            test=True,
                            paired=True,
                            outfile=fcog_outfile,
                            threshold=3)

# taxout = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/rudney3pairs/rudney_3pairs_taxonomy_maxquant_out.tabular'
#
# t = metaquant.metaquant('tax',
#                         file=infile,
#                         sample1_colnames=['Intensity.737NS', 'Intensity.852NS', 'Intensity.867NS'],
#                         sample2_colnames=['Intensity.737WS', 'Intensity.852WS', 'Intensity.867WS'],
#                         test=True,
#                         paired=True,
#                         outfile=taxout,
#                         threshold=3)
#
# # fun tax interaction
# ftout = '/home/caleb/Griffin_lab_work/functional_taxonomic_quant/rudney3pairs/rudney_3pairs_funtax_maxquant_out.tabular'
# ft = metaquant.metaquant('taxfn',
#                          file=infile,
#                          sample1_colnames=['Intensity.737NS', 'Intensity.852NS', 'Intensity.867NS'],
#                          sample2_colnames=['Intensity.737WS', 'Intensity.852WS', 'Intensity.867WS'],
#                          lca_colname="lca",
#                          cog_colname="cog",
#                          test=True,
#                          paired=True,
#                          outfile=ftout,
#                          threshold=2)