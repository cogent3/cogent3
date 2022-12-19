[![PyPI version](https://badge.fury.io/py/cogent3.svg)](https://badge.fury.io/py/cogent3)
[![Downloads](https://pepy.tech/badge/cogent3/month)](https://pepy.tech/project/cogent3)

[![Build Status](https://github.com/cogent3/cogent3/workflows/CI/badge.svg?branch=develop)](https://github.com/cogent3/cogent3/actions?workflow=CI)
[![coverall](https://coveralls.io/repos/github/cogent3/cogent3/badge.svg?branch=develop)](https://coveralls.io/github/cogent3/cogent3?branch=develop)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/cogent3/cogent3.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/cogent3/cogent3/context:python)

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/release/python-370/)
![Using Black Formatting](https://img.shields.io/badge/code%20style-black-000000.svg)

## `cogent3`

`cogent3` is a mature python library for analysis of genomic sequence data. We endeavour to provide a first-class experience within Jupyter notebooks, but the algorithms also support parallel execution on compute systems with 1000's of processors.

## Who is it for?

### Anyone who wants to analyse sequence divergence using robust statistical models

`cogent3` is unique in providing numerous [non-stationary Markov models](http://www.ncbi.nlm.nih.gov/pubmed/25503772) for modelling sequence evolution, [including codon models](https://www.ncbi.nlm.nih.gov/pubmed/28175284). `cogent3` also includes an extensive collection of time-reversible models (again including [novel codon models](https://www.ncbi.nlm.nih.gov/pubmed/19815689)). We have done more than just invent these new methods, we have [established the most robust algorithms](https://www.ncbi.nlm.nih.gov/pubmed/19099591) for their implementation and their [suitability for real data](https://www.ncbi.nlm.nih.gov/pubmed/23935949). Additionally, there are novel signal processing methods focussed on statistical estimation of [integer period signals](https://www.ncbi.nlm.nih.gov/pubmed/21527008).

![nstat](https://cogent3.org/_static/gif/demo-fit-ns.gif)

### Anyone who wants to undertake exploratory genomic data analysis

Beyond our novel methods, `cogent3` provides an extensive suite of capabilities for manipulating and analysing sequence data. You can manipulate sequences by their annotations, e.g.

![annot](https://cogent3.org/_static/gif/demo-annotate.gif)

Plus, you can read standard tabular and biological data formats, perform multiple sequence alignment using any `cogent3` substitution models, phylogenetic reconstruction and tree manipulation, manipulation of tabular data, visualisation of phylogenies and much more.

### Anyone looking for a functional programming style approach to genomic data analysis

Our `cogent3.app` module provides a very different approach to using the library capabilities. Notably, a functional programming style interface lowers the barrier to entry for using `cogent3`'s advanced capabilities. It also supports building pipelines suitable for large-scale analysis. Individuals comfortable with R should find this interface pretty easy to use.

## Installation?

```bash
$ pip install cogent3
```

### Install `extra` -- adds visualisation support

The `extra` group includes python libraries required for visualisation, i.e. [plotly](https://pypi.org/project/plotly/), [kaleido](https://pypi.org/project/kaleido/), [psutil](https://pypi.org/project/psutil/) and [pandas](https://pypi.org/project/pandas/).

```bash
$ pip install "cogent3[extra]"
```

### Install `dev` -- adds `cogent3` development related libraries

The `dev` group includes python libraries required for development of `cogent3`.

```bash
$ pip install "cogent3[dev]"
```

### Install the development version

```bash
$ pip install git+https://github.com/cogent3/cogent3.git@develop#egg=cogent3
```

## Project Information

`cogent3` is released under the BSD-3 license, documentation is at [cogent3.org](https://cogent3.org), while [`cogent3` code is on GitHub](https://github.com/cogent3/cogent3). If you would like to contribute (and we hope you do!), we have created a companion [`c3dev` GitHub](https://github.com/cogent3/c3dev) repo which provides details on how to contribute and some useful tools for doing so.

## Project History

`cogent3` is a descendant of [PyCogent](https://github.com/pycogent/pycogent.github.com). While there is much in common with PyCogent, the amount of change has been substantial, motivating the name change to `cogent3`. This name has been chosen because `cogent` was always the import name (dating back to [PyEvolve in 2004](https://www.ncbi.nlm.nih.gov/pubmed/14706121)) and it's Python 3 only.

Given this history, we are grateful to the multitude of individuals who have made contributions over the years. These individuals are explicitly acknowledged in all the files they contributed to and were co-authors on the original [PyEvolve](https://www.ncbi.nlm.nih.gov/pubmed/14706121) and [PyCogent](https://www.ncbi.nlm.nih.gov/pubmed/17708774) publications.

Compared to PyCogent version 1.9, there has been a massive amount of changes. These include integration of many of the new developments on algorithms and modelling published by the [Huttley lab](https://biology.anu.edu.au/research/groups/huttley-group-bioinformatics-molecular-evolution-genomes) over the last decade. We have also modernised our dependencies. For example, we now use `plotly` for visualisation, `tqdm` for progress bar display, `concurrent.futures` and `mpi4py.futures` for parallel process execution, `nox` and `pytest` for unit testing.
