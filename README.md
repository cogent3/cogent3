
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

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15067121.svg)](https://doi.org/10.5281/zenodo.15067121)


`cogent3` is a mature python library for analysis of genomic sequence data. We endeavour to provide a first-class experience within Jupyter notebooks, but the algorithms also support parallel execution on compute systems with 1000's of processors.

## 📣 Feature Announcements 📣

<details>
  <summary> Migration to new type core objects ‼️ </summary>

We are changing the migration strategy from old type to new type `cogent3` core classes. At present we have old type and new type implementations for sequences, sequence collections, alignments, molecular types, alphabets and genetic codes. Users can select the new classes by specifying `new_type=True` to the functions like `make_aligned_seqs()` or `load_aligned_seqs()`. Alternately, you can do this across all objects by using the `COGENT3_NEW_TYPE` environment variable. We have established that it is not viable to support both old and new types simultaneously. Therefore, **the first release after July 1st 2025 will remove all of the old type classes!** Arguments specific to the old type classes will be deprecated at that point. While this is a major change, we have been using these ourselves consistently and feel confident that the disruption to users should be small. However, we strongly advise all users to migrate now and report any errors. To do this, add the following statement to the top of your scripts.

```python
import os

os.environ["COGENT3_NEW_TYPE"] = "1"
```

</details>

<details>
  <summary> Major advances in our progress towards a fully plugin-based architecture! </summary>

### Cogent3 supports sequence storage plugins 📦🔌🚀

We have implemented the infrastructure to support alternative sequence storage plugins. These provide the backend storage for the new type sequence collections. We have implemented a proof-of-principle plugin [cogent3-h5seqs](https://pypi.org/project/cogent3-h5seqs/) for sequence storage based on the HDF5 format. This allows efficient storage of very large sequence collections (aligned or unaligned). See the readme for that project on how to use it.

### Cogent3 supports sequence format parser and writer plugins 👓✍️🔌

We have implemented the infrastructure to support third-party provision of every bioinformaticians favourite game -- parsing / writing the multitude of sequence file formats.  All builtin format parsers / writers are implemented as plugins. We use third-party versions by default.


### Cogent3 implements plugin hooks 🔌🪝🎉

We have implemented the infrastructure to support hook-style plugins. We have definied a single hook now -- the new type ``Alignment.quick_tree()`` method checks for an external plugin for calculation. The developers of [piqtree](https://pypi.org/project/piqtree) have made the rapid-NJ algorithm available for this hook! Once installed, it is used as `aln.quick_tree(use_hook="piqtree")`.

> **Note**
> For assistance in writing your own plugins, contact us via the [cogent3 discussions page](https://github.com/cogent3/cogent3/discussions).

</details>

<details>
  <summary> Now distributed with sample data! </summary>

  We have added sample data sets for quick testing of different features. Check out `cogent3.available_datasets()` to see the available datasets. You can load one using `cogent3.get_dataset(name)`.

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

## Installation

For most users we recommend

```bash
$ pip install "cogent3[extra]"
```

which installs support for data visualisation and jupyter notebooks.

If you're running on a high-performance computing system we recommend

```bash
$ pip install cogent3
```

which skips the data visualisation and notebook support.

To install the development version directly from GitHub

```bash
$ pip install git+https://github.com/cogent3/cogent3.git@develop#egg=cogent3
```

## Project Information

`cogent3` is released under the BSD-3 license, documentation is at [cogent3.org](https://cogent3.org), while [`cogent3` code is on GitHub](https://github.com/cogent3/cogent3). If you would like to contribute (and we hope you do!), we have created a companion [`c3dev` GitHub](https://github.com/cogent3/c3dev) repo which provides details on how to contribute and some useful tools for doing so.

## Project History

`cogent3` is a descendant of [PyCogent](https://github.com/pycogent/pycogent.github.com). While there is much in common with PyCogent, the amount of change has been substantial, motivating the name change to `cogent3`. This name has been chosen because `cogent` was always the import name (dating back to [PyEvolve in 2004](https://www.ncbi.nlm.nih.gov/pubmed/14706121)) and it's Python 3 only.

Given this history, we are grateful to the multitude of individuals who have made contributions over the years. Many of these contributors were also co-authors on the original [PyEvolve](https://www.ncbi.nlm.nih.gov/pubmed/14706121) and [PyCogent](https://www.ncbi.nlm.nih.gov/pubmed/17708774) publications. Individual contributions can be seen by using "view git blame" on individual lines of code on GitHub, through git log in the terminal, and more recently the changelog.

## Funding

Cogent3 has received funding support from the Australian National University and an [Essential Open Source Software for Science Grant](https://chanzuckerberg.com/eoss/proposals/cogent3-python-apis-for-iq-tree-and-graphbin-via-a-plug-in-architecture/) from the Chan Zuckerberg Initiative.

<p align="center">
  &nbsp;&nbsp;&nbsp;&nbsp;
  <img src="https://webstyle.anu.edu.au/_anu/4/images/logos/2x_anu_logo_small.svg" height="100">
  &nbsp;&nbsp;&nbsp;&nbsp;
  <img src="https://chanzuckerberg.com/wp-content/themes/czi/img/logo.svg" height="110">
</p>
