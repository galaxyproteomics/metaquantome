# metaQuant

Quantitative analysis of the function and taxonomy of microbiomes and their interaction.

# Setup

The newest version of metaquant should be downloaded from this site.
The dependencies are most easily satisfied with conda, and the environment can
be created as follows:

`
conda create -n metaquant python=3.5 pandas ete3 goatools wget numpy statsmodels biopython
`

Note that the bioconda and conda forge channels must be enabled,
as described on [the bioconda website](https://bioconda.github.io/#set-up-channels).

# Tests
To run unittests for the project, run the following from the root directory:

`
python -m unittest discover tests
`


# Roadmap

## High Priority
- visualizations
- documentation
- switch to metagenomics slim
- better arg checking
    - deal with sample that is completely missing values
    - deal with missing values better
    - taxonomy rank checking
    - move ncbi checking to IO
    - check that supplied columns are in the dataframe
    - check that sample info provides a list, or is coerced to a list
    - don't return above phylum
    - don't return BP, MF, or CC
    - raise error when all rows are filtered out
    - strip any newlines in samps file
- configure Travis CI on Github
- add option for specific rank in TF
- move threshold to calculating mean values as well
- make COG more similar to other databases
    - make class
    - add is_in_db method

## Lower Priority
- benchmarking and optimization
- use flake8 for codestyle
- refactoring
    - move every analysis function to one to cut down on duplication?
    - could implement databases as instantiations of abstract classes. Any advantages?

## Done
- unified database structure and 'adding up'
    - unify this with classes
    - work on EC in particular
    - how do we implement for COG cats?
- add handling if description is not found in database (EC)

