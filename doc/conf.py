import os


extensions = [
    "numpydoc",
    "sphinx.ext.todo",
    "sphinx.ext.doctest",
    "nbsphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinxcontrib.bibtex",
]

# todo_include_todos=True # to expose the TODOs, uncomment this line

# Add any paths that contain templates here, relative to this directory.
templates_path = ["templates"]

# The suffix of source filenames.
source_suffix = ".rst", ".ipynb"

# ignore the cookbook/union_dict.rst file as it's specifically included
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints", "cookbook/union_dict.rst"]

# The encoding of source files.
# source_encoding = 'utf-8'

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "Cogent3"
copyright = "2019, Cogent3"

version = ""

release = "2019.10.17a"

# exclude_trees = ["_build"]

show_authors = True

pygments_style = "sphinx"

# -- Options for HTML output ---------------------------------------------------

html_theme = "alabaster"

html_static_path = ["_static"]

htmlhelp_basename = "Cogent3doc"


# -- Options for LaTeX output --------------------------------------------------
latex_documents = [
    ("index", "Cogent3.tex", "Cogent3 Documentation", "Cogent3 Team", "manual")
]


def setup(app):
    app.add_javascript("require.min.js")
    app.add_javascript("plotly-latest.min.js")
