import tax
import matplotlib.pyplot as plt

trey = tax.TreeOfLife('test/unipept_results_metaproteomics.tabular', ['intensity'])

abund = trey.rel_abundance_rank('phylum').plot(kind = "bar", legend = None)
plt.subplots_adjust(bottom = 0.45)
plt.show()

