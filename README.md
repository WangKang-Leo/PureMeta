## PureMeta

## Introduction

[PureMeta](https://ecotyper.stanford.edu/) is an integrated bioinfomatic workflow to qunantify the metabolic states of tumor cells from bulk gene expression data.

We have already defined cell states and ecotypes across **carcinomas** ([Luca/Steen et al., Cell 2021](https://doi.org/10.1016/j.cell.2021.09.014)) and in **diffuse large B cell lymphoma (DLBCL)** ([Steen/Luca et al., Cancer Cell 2021](https://doi.org/10.1016/j.ccell.2021.08.011)). The current version of EcoTyper allows users to recover the cell states and ecotypes for these two tumor categories in their own data. Additionally, it allows users to discover and recover cell states and ecotypes in their system of interest, including **directly** from scRNA-seq data (see [Tutorial 5](#tutorial-5-de-novo-discovery-of-cell-states-and-ecotypes-in-scrna-seq-data)). Below we illustrate each of these functionalities. 

## Citation

If EcoTyper software, data, and/or website are used in your publication, please cite the following paper(s): 

* [Luca/Steen et al., Cell 2021](https://doi.org/10.1016/j.cell.2021.09.014) (detailed description of EcoTyper and application to carcinomas).
* [Steen/Luca et al., Cancer Cell 2021](https://doi.org/10.1016/j.ccell.2021.08.011) (application of EcoTyper to lymphoma).

## Setup

The latest version of EcoTyper source code can be found on [EcoTyper GitHub repository](https://github.com/digitalcytometry/ecotyper) and [Ecotyper website](https://ecotyper.stanford.edu/). To set up EcoTyper, please download this folder locally:

```{bash, eval = F}
git clone https://github.com/digitalcytometry/ecotyper
cd ecotyper
```

or:

```{bash, eval = F}
wget https://github.com/digitalcytometry/ecotyper/archive/refs/heads/master.zip
unzip master.zip
cd ecotyper-master
```
