
<p align="left">
  <img src="https://raw.githubusercontent.com/cogent3/cogent3.github.io/e72df8c155c100f502b6a7009347d1821ab3adef/doc/_static/c3-logo.svg" width="300">
</p>

[![PyPI version](https://badge.fury.io/py/cogent3.svg)](https://badge.fury.io/py/cogent3)
[![PyPI Downloads](https://pepy.tech/badge/cogent3/month)](https://pepy.tech/project/cogent3)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cogent3/badges/downloads.svg)](https://anaconda.org/bioconda/cogent3)

[![Build Status](https://github.com/cogent3/cogent3/workflows/CI/badge.svg?branch=develop)](https://github.com/cogent3/cogent3/actions?workflow=CI)
[![coverall](https://coveralls.io/repos/github/cogent3/cogent3/badge.svg?branch=develop)](https://coveralls.io/github/cogent3/cogent3?branch=develop)

![PyPI - Python Version](https://img.shields.io/pypi/pyversions/cogent3)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![CodeQL](https://github.com/cogent3/cogent3/actions/workflows/codeql.yml/badge.svg)](https://github.com/cogent3/cogent3/actions/workflows/codeql.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e80e3441de59449bb1a4d8ad1fdea4fa)](https://app.codacy.com/gh/cogent3/cogent3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15067121.svg)](https://doi.org/10.5281/zenodo.15067121)

[![CZI's Essential Open Source Software for Science](https://img.shields.io/badge/funded%20by-EOSS-FF414B)](https://czi.co/EOSS)

`cogent3` is a mature python library for analysis of genomic sequence data. We endeavour to provide a first-class experience within Jupyter notebooks, but the algorithms also support parallel execution on compute systems with 1000's of processors. A plugin system enables custom extensions to expand the library capabilities.

## 📣 Features & Announcements 📣

<details>
  <summary>App foundations and more have a new home -- scinexus!</summary>

The `cogent3` apps and plugin architecture remain the same, but we have extracted the app infrastructure into a new project called [scinexus](https://pypi.org/project/scinexus/). The objective is to make this powerful framework available beyond cogent3. We also migrated the data stores and many generally useful utilities. e.g. `open_`, deserialisation of objects, and progress bars. Check out the documentation of scinexus to see what it provides.

</details>

<details>
  <summary>Tracking citations of apps</summary>

We have built a mechanism for defining and tracking citations for apps, making it easier for users to correctly acknowledge app developer efforts. See the announcement on the home page for more details.

</details>

<details>
  <summary>Drawing annotations from an annotation database</summary>

The new function `cogent3.draw_annotations()` can take an annotation database instance and returns a Drawable. This is useful for exploration of genome features without needing to load sequences.

</details>

<details>
  <summary>The cogent3 code-sharing site 📚🔌</summary>

Share your [cogent3 ecosystem code solutions](https://github.com/cogent3/c3codeshare) for others to benefit from your awesomeness 😎. Or, just read / use the contributions from others.

</details>

<details>
  <summary> We now support third-party plugins for annotation databases 📚🔌 </summary>

If you want to contribute a third-party backend, [get in touch](https://github.com/cogent3/cogent3/discussions)!

</details>

<details>
  <summary> diverse-seq has been rewritten in rust 🚀! </summary>

The sequence sampling tool [diverse-seq](diverse-seq.readthedocs.io), which provides multiple apps for sampling representative sequences, just got faster! The performance critical code has been rewritten in Rust. Give it a try 😀.

</details>

<details>
  <summary> A new rust-based plugin for k-mer counting </summary>

We recently added a new `count_kmers()` method to the `SequenceCollection` and `Sequence` classes. Then, the developers of [Pykmertools](https://github.com/anuradhawick/kmertools) (with a bit of help from us) have released a `cogent3-pykmertools` app which makes their rust-based python module for counting k-mers available as `seqs.count_kmers(k=k, use_hook="cogent3_pykmertools")`. Install it with `pip install cogent3-pykmertools` and give it a try. And add a star to the [Pykmertools](https://github.com/anuradhawick/kmertools) repo!

</details>

<details>
  <summary> Major advances in our progress towards a fully plugin-based architecture! </summary>

### Cogent3 supports sequence storage plugins 📦🔌🚀

We have implemented the infrastructure to support alternative sequence storage plugins. These provide the backend storage for the new type sequence collections. We have implemented a proof-of-principle plugin [cogent3-h5seqs](https://pypi.org/project/cogent3-h5seqs/) for sequence storage based on the HDF5 format. This allows efficient storage of very large sequence collections (aligned or unaligned). See the readme for that project on how to use it.

### Cogent3 supports sequence format parser and writer plugins 👓✍️🔌

We have implemented the infrastructure to support third-party provision of every bioinformaticians favourite game -- parsing / writing the multitude of sequence file formats.  All builtin format parsers / writers are implemented as plugins. We use third-party versions by default.

### Cogent3 implements plugin hooks 🔌🪝🎉

We have implemented the infrastructure to support hook-style plugins. We have defined a single hook now -- the new type `Alignment.quick_tree()` method checks for an external plugin for calculation. The developers of [piqtree](https://pypi.org/project/piqtree) have made the rapid-NJ algorithm available for this hook! Once installed, it is used as `aln.quick_tree(use_hook="piqtree")`.

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

which installs support for data visualisation (which requires [Plotly](https://pypi.org/project/plotly/)) and extensions for jupyter notebooks.

> **Note:** In order to write Plotly figures to static image files you will need to [install Chrome](https://plotly.com/python/static-image-export/#install-dependencies).

### Minimal installation

If you're running on a high-performance computing system we recommend

```bash
$ pip install cogent3
```

which skips the data visualisation and notebook support.

### Install with developer tools

Everything we use for `cogent3` development.

```bash
$ pip install "cogent3[dev]"
```

> **Note:** Installs all dependencies that can be installed using `pip`.

### Installing the development version

```bash
$ pip install git+https://github.com/cogent3/cogent3.git@develop#egg=cogent3
```

> **Warning:** The interface can change without warning on the development branch.

### Installing using conda / mamba

Activate your [conda](https://docs.conda.io/en/latest/miniconda.html) environment, then

```bash
(myenv) $ conda install bioconda::cogent3
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
