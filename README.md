[![Build Status](https://travis-ci.org/galaxyproteomics/metaquantome.svg?branch=master)](https://travis-ci.org/galaxyproteomics/metaquantome)

[![Coverage Status](https://coveralls.io/repos/github/galaxyproteomics/metaquantome/badge.svg?branch=master)](https://coveralls.io/github/galaxyproteomics/metaquantome?branch=master)

# metaQuantome


Quantitative analysis on the function and taxonomy of microbiomes and their interaction.

## Table of Contents
- [Setup](#setup-)
    - [BioConda](#bioconda-)
    - [Local Development](#local-development-)
- [Usage](#usage-)
- [Tests](#tests-)

## Setup (#table-of-contents)

### BioConda

The easiest way to install metaQuantome with all the dependencies is by using Bioconda (provided you are on Mac or Linux, which are the only systems supported by Bioconda).

First, install the conda package manager, by downloading either Anaconda or Miniconda (see https://docs.anaconda.com/anaconda/install/). Then, the following commands will set up the necessary channels for Bioconda and metaQuantome. If needed, additional information on Bioconda is available at https://bioconda.github.io/.

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Then, run the following command to set up an environment named mqome, which will have metaQuantome (version 2.0.1) and all dependencies in it:

```
conda create -n mqome metaquantome=2.0.1
```

If the following prompt is seen at the command line, type `y`:
```
Proceeed ([y]/n)?
```

The create command only needs to be run once.

Finally, we can activate the environment using the following command, which will make the `metaquantome` command available.
```
conda activate mqome
```

### Local Development

To begin, ensure that you have Python3 installed. To check, issue this command to verify your python version:
```
python --version
```

If Python3 is not installed, please download it from [here](https://www.python.org/downloads/).

Pip is the package installer for Python. It comes pre-packaged with Python. 

To install your local version of metaQuantome, run the following in the root directory:
```
pip install .
```

## Usage (#table-of-contents)

In-depth tutorials can be found [here](https://galaxyproteomics.github.io/metaquantome_mcp_analysis/) for the following:
1. [metaQuantome command-line interface](https://galaxyproteomics.github.io/metaquantome_mcp_analysis/cli_tutorial/cli_tutorial.html)
2. [metaQuantome Galaxy tool](https://galaxyproteomics.github.io/metaquantome_mcp_analysis/galaxy_tutorial/galaxy_tutorial.html)

## Tests (#table-of-contents)


For testing, you will also need to set up a BioConda environment as specified [above](#bioconda-).

Before running unittests for the project, you will need to install some databases as follows:
```
metaquantome db ncbi --dir ./metaquantome/data/test/
```
*Note: This step will take a few minutes to download*

Once set up, run the following from the root directory:

```sh
python -m unittest discover tests
```
*Note: This step will take a few minutes to run*
