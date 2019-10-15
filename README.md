[![Build Status](https://dev.azure.com/GavinHuttley/cogent3/_apis/build/status/cogent3.cogent3?branchName=master)](https://dev.azure.com/GavinHuttley/cogent3/_build/latest?definitionId=1&branchName=master)
[![codecov](https://codecov.io/gh/cogent3/cogent3/branch/master/graph/badge.svg)](https://codecov.io/gh/cogent3/cogent3)
![Using Black Formatting](https://img.shields.io/badge/code%20style-black-000000.svg)

## `cogent3`

`cogent3` is a python library for genomic data analysis, specifically comparative genomics data.

## Installation?

We will be up on PyPI as soon as the documentation is written, until then it's [download a tarball](https://github.com/cogent3/cogent3) and pip install from that, or clone.

## Project Information

`cogent3` is released under the BSD-3 license, documentation for [`cogent3` is on readthedocs](https://cogent3.readthedocs.io/en/latest/), while [`cogent3` code is GitHub](https://github.com/cogent3/cogent3). If you would like to contribute (and we hope you do!), we have created a companion [`c3dev` GitHub](https://github.com/cogent3/c3dev) repo which provides details on how to contribute and some useful tools for doing so.

## Project History

`cogent3` is a descendant of [PyCogent](https://github.com/pycogent/pycogent.github.com). While there is much in common with PyCogent, the amount of change has been substantial, motivating a new name `cogent3`. This name has been chosen because `cogent` was always the import name (dating back to 2004!) and it's Python 3 only.

Compared to PyCogent version 1.9, there have been a massive amount of changes. These include integration of many of the new developments on modelling of non-stationary processes published by the Huttley lab over the last decade. We have also modernised some convenience capabilities. For example, we now use `tqdm` for progress bar display and `concurrent.futures` and `mpi4py.futures` for parallel process execution.

We have implemented a `cogent3.app` module which contains a very different approach to using the capabilities. Notably, a functional programming style interface lowers the barrier to entry for using `cogent3`'s advanced capabilities. Documentation is coming...

Importantly, the app interface should be considered as alpha level code.
