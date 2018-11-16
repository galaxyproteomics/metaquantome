library(ggplot2)

# source the read data function
source('read_result.R')


ggplot(ns) +
    geom_bar(aes(x = reorder(taxon_name, -NS_mean), y = 2^NS_mean), stat = "identity") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
    labs(x = "Taxon", y = "Total Peptide Intensity",
         title = "Most Abundant Genera in NS")
