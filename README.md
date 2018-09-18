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
- better arg checking
- configure Travis CI on Github
- add option for specific rank in TF

## Lower Priority
- use flake8 for codestyle

## Done
- unified database structure and 'adding up'
    - unify this with classes
    - work on EC in particular
    - how do we implement for COG cats?
- add handling if description is not found in database (EC)

