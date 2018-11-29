[![Build Status](https://travis-ci.org/galaxyproteomics/metaquantome.svg?branch=master)](https://travis-ci.org/galaxyproteomics/metaquantome)

[![Coverage Status](https://coveralls.io/repos/github/galaxyproteomics/metaquantome/badge.svg?branch=master)](https://coveralls.io/github/galaxyproteomics/metaquantome?branch=master)

# metaQuantome


Quantitative analysis of the function and taxonomy of microbiomes and their interaction.

# Setup

The newest version of metaquantome should be downloaded from this site.
The dependencies are most easily satisfied with conda, and the environment can
be created and activated as follows:

```sh
conda env update --file dev_environment.yml
source activate metaquantome
```

Note that the bioconda and conda forge channels must be enabled,
as described on [the bioconda website](https://bioconda.github.io/#set-up-channels).

# Tests
To run unittests for the project, run the following from the root directory:

```sh
python -m unittest discover tests
```

# Basic usage

We can use example files within the `examples` directory to demonstrate the use of the command line.

```sh
python3 cli.py expand -m f \
    --pep_colname peptide --outfile test.tmp -i example/int.tab \
    --func_file example/func.tab  --func_colname go --ontology go \
    --samps '{"s1": ["int1", "int2", "int3"], "s2": ["int4", "int5", "int6"]}'
```

# Roadmap

## Critical changes before Bioconda release
- change modes to 't', 'f', 'ft'
- add peptide column name for each type of file


## High Priority
- visualizations
- documentation
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


## Lower Priority
- benchmarking and optimization
- use flake8 for codestyle

## Done
- unified database structure and 'adding up'
    - unify this with classes
    - work on EC in particular
    - how do we implement for COG cats?
- add handling if description is not found in database (EC)
- switch to metagenomics slim
- add option for specific rank in TF
- move threshold to calculating mean values as well
