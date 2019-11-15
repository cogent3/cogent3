[![Build Status](https://dev.azure.com/GavinHuttley/cogent3/_apis/build/status/cogent3.cogent3?branchName=master)](https://dev.azure.com/GavinHuttley/cogent3/_build/latest?definitionId=1&branchName=master)
[![codecov](https://codecov.io/gh/cogent3/cogent3/branch/master/graph/badge.svg)](https://codecov.io/gh/cogent3/cogent3)
![Using Black Formatting](https://img.shields.io/badge/code%20style-black-000000.svg)

## `cogent3`

`cogent3` is a mature python library for analysis of genomic sequence data. We endeavour to provide a first-class experience within Jupyter notebooks, but the algorithms also support parallel execution on compute systems with 1000's of processors.

## Who is it for?

### Anyone who wants to analyse sequence divergence using robust statistical models

`cogent3` is unique in providing numerous [non-stationary Markov models](http://www.ncbi.nlm.nih.gov/pubmed/25503772) for modelling sequence evolution, [including codon models](https://www.ncbi.nlm.nih.gov/pubmed/28175284). `cogent3` also includes an extensive collection of time-reversible models (again including [novel codon models](https://www.ncbi.nlm.nih.gov/pubmed/19815689)). We have done more than just invent these new methods, we have [established the most robust algorithms](https://www.ncbi.nlm.nih.gov/pubmed/19099591) for their implementation and their [suitability for real data](https://www.ncbi.nlm.nih.gov/pubmed/23935949). Additionally, there are novel signal processing methods focussed on statistical estimation of [integer period signals](https://www.ncbi.nlm.nih.gov/pubmed/21527008).

### Anyone who wants to undertake exploratory genomic data analysis

Beyond our novel methods, `cogent3` provides an extensive suite of capabilities for manipulating and analysing sequence data. For instance, the ability to read standard biological data formats, manipulate sequences by their annotations, to perform multiple sequence alignment using any of our substitution models, phylogenetic reconstruction and tree manipulation, manipulation of tabular data, visualisation of phylogenies and much more.

### Anyone looking for a functional programming style approach to genomic data analysis

Our `cogent3.app` module provides a very different approach to using the library capabilities. Notably, a functional programming style interface lowers the barrier to entry for using `cogent3`'s advanced capabilities. It also supports building pipelines suitable for large-scale analysis. Individuals comfortable with R should find this interface pretty easy to use.

## Installation?

```bash
$ pip install cogent3
```

### Installing the development version

```bash
$ pip install git+https://github.com/cogent3/cogent3.git@master#egg=cogent3
```

## Project Information

`cogent3` is released under the BSD-3 license, documentation for [`cogent3` is on readthedocs](https://cogent3.readthedocs.io/en/latest/), while [`cogent3` code is on GitHub](https://github.com/cogent3/cogent3). If you would like to contribute (and we hope you do!), we have created a companion [`c3dev` GitHub](https://github.com/cogent3/c3dev) repo which provides details on how to contribute and some useful tools for doing so.

## Project History

`cogent3` is a descendant of [PyCogent](https://github.com/pycogent/pycogent.github.com). While there is much in common with PyCogent, the amount of change has been substantial, motivating a new name `cogent3`. This name has been chosen because `cogent` was always the import name (dating back to [PyEvolve in 2004](https://www.ncbi.nlm.nih.gov/pubmed/14706121)) and it's Python 3 only.

Given this history, we are grateful to the multitude of individuals who have made contributions over the years. These individuals are explicitly acknowledged in all the files they contributed to and were co-authors on the original [PyEvolve](https://www.ncbi.nlm.nih.gov/pubmed/14706121) and [PyCogent](https://www.ncbi.nlm.nih.gov/pubmed/17708774) publications.

Compared to PyCogent version 1.9, there has been a massive amount of changes. These include integration of many of the new developments on algorithms and modelling published by the Huttley lab over the last decade. We have also modernised our dependencies. For example, we now use `plotly` for visualisation, `tqdm` for progress bar display, `concurrent.futures` and `mpi4py.futures` for parallel process execution, `tox` and `pytest` for unit testing.
