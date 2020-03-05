import os
import sys


sys.path.append("../src")

# Allow autosummary to generate stub files
autosummary_generate = True

# Prevent numpydoc from requiring stub files for methods
numpydoc_class_members_toctree = False

extensions = [
    "numpydoc",
    "sphinx.ext.todo",
    "sphinx.ext.doctest",
    "nbsphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinxcontrib.bibtex",
    "sphinx_gallery.gen_gallery",
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
project = "cogent3"
copyright = "2020, cogent3"

release = "2020.2.7a"

version = release

# exclude_trees = ["_build"]

show_authors = True

pygments_style = "sphinx"

# -- Options for HTML output ---------------------------------------------------

html_theme = "alabaster"

sidebar_collapse = True

html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
        # 'donate.html',
    ]
}


html_theme_options = {
    "fixed_sidebar": True,
}

html_static_path = ["_static"]

htmlhelp_basename = "cogent3doc"


# -- Options for LaTeX output --------------------------------------------------
latex_documents = [
    ("index", "cogent3.tex", "cogent3 Documentation", "cogent3 Team", "manual")
]

nbsphinx_requirejs_path = "require.js"

from plotly.io._sg_scraper import plotly_sg_scraper
def plotly_sg_scraper_nb(*args, **kwargs):
    try:
        return plotly_sg_scraper(*args, **kwargs)
    except Exception:
        return ""

image_scrapers = ("matplotlib", plotly_sg_scraper_nb,)

from sphinx_gallery.sorting import ExplicitOrder

sphinx_gallery_conf = {
     "doc_module": ("plotly",),
     "examples_dirs": ["draw_examples"],
     "subsection_order": ExplicitOrder(["draw_examples/alignments",
                                       "draw_examples/trees"]),
     "gallery_dirs": ["draw"],
     "backreferences_dir": "api/draw",
     "reference_url": {"plotly": None,
      },
     "image_scrapers": image_scrapers,
}
