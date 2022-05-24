import datetime
import os
import shutil
import pathlib
import sys

from glob import glob

import sphinx_bootstrap_theme

sys.path.append("../src")


def make_nbsphinx_thumbnails():
    """returns dict of {path: '_images/{path.stem}'"""
    gallery = sorted(
        p
        for p in pathlib.Path("draw").glob("**/*.rst")
        if p.stem not in ("index", "README")
    )

    return {str(n).split(".")[0]: f"_images/{n.stem}.png" for n in gallery}


# Allow autosummary to generate stub files
autosummary_generate = True
add_module_names = False  # don't include module path to module/func names

# Prevent numpydoc from requiring stub files for methods
numpydoc_class_members_toctree = False

extensions = [
    "jupyter_sphinx",
    "nbsphinx",
    "numpydoc",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.githubpages",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx_gallery.load_style",
    "sphinx_panels",
    # "sphinxcontrib.spelling",
]

# todo_include_todos=True # to expose the TODOs, uncomment this line

# Add any paths that contain templates here, relative to this directory.
templates_path = ["templates"]

# The suffix of source filenames.
source_suffix = ".rst", ".ipynb"

# ignoring the cookbook/union_dict.rst file as it's specifically included
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "cookbook/union_dict.rst",
    "cookbook/loading_tabular",
    "COGENT3_LICENSE.rst",
    "*tmp*",
]

# The encoding of source files.

# The master toctree document.
master_doc = "index"

# General information about the project.
today = datetime.date.today()
year = today.strftime("%Y")
project = "cogent3"
copyright = f"2020-{year}, cogent3 Team"

release = "2022.5.25a1"

version = ""

# exclude_trees = ["_build"]

show_authors = True

pygments_style = "sphinx"

# -- Options for HTML output ---------------------------------------------------
html_theme = "bootstrap"
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()


html_theme_options = {
    "navbar_title": "Docs",
    "navbar_site_name": "Sections",
    "navbar_links": [
        ("Install", "install"),
        ("Gallery", "draw/index.html", True),
    ],
    "navbar_class": "navbar navbar-inverse",
    "navbar_fixed_top": "true",
    "source_link_position": "skipped",
    "bootswatch_theme": "cerulean",
    "bootstrap_version": "3",
}

html_static_path = ["_static"]

htmlhelp_basename = "cogent3doc"

# -- Options for Gallery

nbsphinx_requirejs_path = "require.js"
nbsphinx_thumbnails = make_nbsphinx_thumbnails()


# -- Options for LaTeX output --------------------------------------------------
latex_documents = [
    ("index", "cogent3.tex", "cogent3 Documentation", "cogent3 Team", "manual")
]


def setup(app):
    app.add_js_file("plotly-latest.min.js")
