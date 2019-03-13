[![Build Status](https://travis-ci.org/galaxyproteomics/metaquantome.svg?branch=master)](https://travis-ci.org/galaxyproteomics/metaquantome)

[![Coverage Status](https://coveralls.io/repos/github/galaxyproteomics/metaquantome/badge.svg?branch=master)](https://coveralls.io/github/galaxyproteomics/metaquantome?branch=master)

# metaQuantome


Quantitative analysis of the function and taxonomy of microbiomes and their interaction.


# Developers

## Setup

The newest version of metaquantome should be downloaded from this site.
The dependencies are most easily satisfied with conda, and the environment can
be created and activated as follows:

```sh
conda env update --file dev_environment.yml
source activate metaquantome
```

Note that the bioconda and conda forge channels must be enabled,
as described on [the bioconda website](https://bioconda.github.io/#set-up-channels).

## Tests
To run unittests for the project, run the following from the root directory:

```sh
python -m unittest discover tests
```
