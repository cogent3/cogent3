
<p align="left">
  <img src="https://raw.githubusercontent.com/cogent3/cogent3.github.io/e72df8c155c100f502b6a7009347d1821ab3adef/doc/_static/c3-logo.svg" width="300">
</p>

[![PyPI version](https://badge.fury.io/py/cogent3.svg)](https://badge.fury.io/py/cogent3)
[![Downloads](https://pepy.tech/badge/cogent3/month)](https://pepy.tech/project/cogent3)

[![Build Status](https://github.com/cogent3/cogent3/workflows/CI/badge.svg?branch=develop)](https://github.com/cogent3/cogent3/actions?workflow=CI)
[![coverall](https://coveralls.io/repos/github/cogent3/cogent3/badge.svg?branch=develop)](https://coveralls.io/github/cogent3/cogent3?branch=develop)

![PyPI - Python Version](https://img.shields.io/pypi/pyversions/cogent3)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![CodeQL](https://github.com/cogent3/cogent3/actions/workflows/codeql.yml/badge.svg)](https://github.com/cogent3/cogent3/actions/workflows/codeql.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e80e3441de59449bb1a4d8ad1fdea4fa)](https://app.codacy.com/gh/cogent3/cogent3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

`cogent3` is a mature python library for analysis of genomic sequence data. We endeavour to provide a first-class experience within Jupyter notebooks, but the algorithms also support parallel execution on compute systems with 1000's of processors.

## 📣 Feature Announcements 📣

<details>
    <summary> Faster sequence coevolution measures 🚀 </summary>

We have completely rewritten all the Mutual Information based coevolution statistic calculators. Single CPU performance is orders of magnitude faster than the old implementation and we now also support parallel execution. The existing `Alignment.coevolution()` method uses these so you don't need to do anything different to use the new algorithms.

</details>

<details>
    <summary> Faster sequence format parsers 💨 </summary>

We have faster implementations of the parsers for Fasta and GenBank sequence formats. These are used by our standard loading mechanisms. If you just want to get the contents of files in those formats as standard Python types, use `cogent3.parser.fasta.iter_fasta_records()` or `cogent3.parser.genbank.iter_genbank_records()`.

</details>

<details>
  <summary> Supporting third-party apps as plugins 🔌 </summary>

Cogent3 now provides support for plugins! Third-party developers can deploy their code as cogent3 apps with just a few lines. See the [demo project](https://github.com/cogent3/app_template).

Post any questions you have in [cogent3 discussions](https://github.com/cogent3/cogent3/discussions).

</details>

<details>
  <summary> The developers of Cogent3 and IQ-TREE2 announce piqtree2 🎉 </summary>

Speaking of plugins, our first major third-party plugin is [piqtree2](https://pypi.org/project/piqtree2). Try it out and [give us feedback](https://github.com/iqtree/piqtree2/discussions).

</details>

<details>
  <summary> New core data types improve efficiency and flexibility </summary>

The cogent3 development team 👾 are hard at work modernising the core internals 💪🛠.

In this release the `Sequence`, `SequenceCollection`, `MolType`, `GeneticCode`, and alphabet classes have all been rewritten from scratch with an eye to simplifying the code while improving its flexibility and performance. (We're working on alignments for the next release.)

The "new-style" objects enhance performance by supporting the access of the underlying data in various formats (i.e. numpy arrays, bytes or strings). You can create "new-style" objects by setting the `new_type=True` argument in top-level functions (`make_seq`, `load_seq`, `make_unaligned_seqs`, `get_moltype`, `get_code`). These are not yet the default and are not fully integrated into the existing code. They can also differ in their API relative to the classes they replace. 

We encourage experimentation in cases where integration with old objects is NOT required and [look forward to any feedback](https://github.com/cogent3/cogent3/discussions)!

</details>

## Who is it for?

### Anyone who wants to analyse sequence divergence using robust statistical models

`cogent3` is unique in providing numerous [non-stationary Markov models](http://www.ncbi.nlm.nih.gov/pubmed/25503772) for modelling sequence evolution, [including codon models](https://www.ncbi.nlm.nih.gov/pubmed/28175284). `cogent3` also includes an extensive collection of time-reversible models (again including [novel codon models](https://www.ncbi.nlm.nih.gov/pubmed/19815689)). We have done more than just invent these new methods, we have [established the most robust algorithms](https://www.ncbi.nlm.nih.gov/pubmed/19099591) for their implementation and their [suitability for real data](https://www.ncbi.nlm.nih.gov/pubmed/23935949). Additionally, there are novel signal processing methods focussed on statistical estimation of [integer period signals](https://www.ncbi.nlm.nih.gov/pubmed/21527008).

<details>
  <summary> 🎬 Demo non-reversible substitution model </summary>
    <video src="https://user-images.githubusercontent.com/3102996/253845402-f511af2c-c2e2-48bc-8f6e-f9b0f05697e9.mp4" controls="controls" style="max-height:640px">
    </video>
</details>

### Anyone who wants to undertake exploratory genomic data analysis

Beyond our novel methods, `cogent3` provides an extensive suite of capabilities for manipulating and analysing sequence data. You can manipulate sequences by their annotations, e.g.

<details>
  <summary> 🎬 Demo sequences with annotations </summary>
    <video src="https://user-images.githubusercontent.com/3102996/253847297-2611cda8-e078-4b86-a269-43fbf6ced14c.mp4" controls="controls" style="max-height:640px">
    </video>
</details>

Plus, you can read standard tabular and biological data formats, perform multiple sequence alignment using any `cogent3` substitution models, phylogenetic reconstruction and tree manipulation, manipulation of tabular data, visualisation of phylogenies and much more.

### Beginner friendly approach to genomic data analysis

Our `cogent3.app` module provides a very different approach to using the library capabilities. Expertise in structural programming concepts is not essential!

<details>
  <summary> 🎬 Demo friendly coding </summary>
    <video src="https://user-images.githubusercontent.com/3102996/253849168-a821de1a-1aad-4761-970f-e365f6b3b1cd.mp4" controls="controls" style="max-height:640px">
    </video>
</details>

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

Given this history, we are grateful to the multitude of individuals who have made contributions over the years. Many of these contributors were also co-authors on the original [PyEvolve](https://www.ncbi.nlm.nih.gov/pubmed/14706121) and [PyCogent](https://www.ncbi.nlm.nih.gov/pubmed/17708774) publications. Individual contributions can be seen by using "view git blame" on individual lines of code on GitHub , through git log in the terminal, and more recently the changelog.

Compared to PyCogent version 1.9, there has been a massive amount of changes. These include integration of many of the new developments on algorithms and modelling published by the [Huttley lab](https://biology.anu.edu.au/research/groups/huttley-group-bioinformatics-molecular-evolution-genomes) over the last decade. We have also modernised our dependencies. For example, we now use `plotly` for visualisation, `tqdm` for progress bar display, `concurrent.futures` and `mpi4py.futures` for parallel process execution, `nox` and `pytest` for unit testing.

## Funding

Cogent3 has received funding support from the Australian National University and an [Essential Open Source Software for Science Grant](https://chanzuckerberg.com/eoss/proposals/cogent3-python-apis-for-iq-tree-and-graphbin-via-a-plug-in-architecture/) from the Chan Zuckerberg Initiative.


<p align="center">
  &nbsp;&nbsp;&nbsp;&nbsp;
  <img src="https://webstyle.anu.edu.au/_anu/4/images/logos/2x_anu_logo_small.svg" height="100">
  &nbsp;&nbsp;&nbsp;&nbsp;
  <img src="https://chanzuckerberg.com/wp-content/themes/czi/img/logo.svg" height="110">
</p>
