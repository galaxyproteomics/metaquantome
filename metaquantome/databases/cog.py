def take_first_cog(df, cog_name):
    """
    In a list of COGs, assume the first is the most relevant
    todo: could change this to normalize the dataframe, as in other ontologies

    :param df: DataFrame
    :param cog_name: column with COGs
    :return: dataframe with one COG per row.
    """
    # take first cog
    df[cog_name] = df[cog_name].str.split(',').str[0]
    return df


# from ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab
# mapping letters to description
cogCat = {'J': 'Translation, ribosomal structure and biogenesis',
          'A': 'RNA processing and modification',
          'K': 'Transcription',
          'L': 'Replication, recombination and repair',
          'B': 'Chromatin structure and dynamics',
          'D': 'Cell cycle control, cell division, chromosome partitioning',
          'Y': 'Nuclear structure',
          'V': 'Defense mechanisms',
          'T': 'Signal transduction mechanisms',
          'M': 'Cell wall/membrane/envelope biogenesis',
          'N': 'Cell motility',
          'Z': 'Cytoskeleton',
          'W': 'Extracellular structures',
          'U': 'Intracellular trafficking, secretion, and vesicular transport',
          'O': 'Posttranslational modification, protein turnover, chaperones',
          'X': 'Mobilome: prophages, transposons',
          'C': 'Energy production and conversion',
          'G': 'Carbohydrate transport and metabolism',
          'E': 'Amino acid transport and metabolism',
          'F': 'Nucleotide transport and metabolism',
          'H': 'Coenzyme transport and metabolism',
          'I': 'Lipid transport and metabolism',
          'P': 'Inorganic ion transport and metabolism',
          'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
          'R': 'General function prediction only',
          'S': 'Function unknown',
          'unknown': 'unassigned function'}
