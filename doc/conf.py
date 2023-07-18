import datetime
import pathlib


def make_nbsphinx_thumbnails():
    """returns dict of {path: '_images/{path.stem}'"""
    gallery = [
        p for p in pathlib.Path("doc/draw").glob("**/*.rst") if p.stem != "README"
    ]

    return {str(n).split(".")[0]: f"_images/{n.stem}.png" for n in gallery}


# sphinx_navtree
today = datetime.date.today()
year = today.strftime("%Y")
project = "cogent3"
copyright = f"2020-{year}, cogent3 Team"
author = "Gavin Huttley"

# The full version, including alpha/beta/rc tags
# Use calendar versioning
release = today.strftime("%Y.%m.%d")


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

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
    "sphinxcontrib.bibtex",
    "sphinx_panels",
]

# stop sphinx-panels from including css
# panels_add_bootstrap_css = False


# Allow autosummary to generate stub files
autosummary_generate = True
add_module_names = False  # don't include module path to module/func names
# Prevent numpydoc from requiring stub files for methods
numpydoc_class_members_toctree = False

bibtex_bibfiles = ["cogent3.bib"]

templates_path = ["doc/templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The master toctree document.
master_doc = "index"
show_authors = True
pygments_style = "sphinx"

todo_include_todos = False
todo_emit_warnings = True
htmlhelp_basename = "cogent3doc"


# ignoring the cookbook/union_dict.rst file as it's specifically included
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "cookbook/union_dict",
    "cookbook/loading_tabular",
    "COGENT3_LICENSE",
    "*tmp*",
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

sidebar_collapse = False

html_theme_options = {
    "navigation_depth": 6,
    "show_toc_level": 4,
    "show_nav_level": 0,
    "github_url": "https://github.com/cogent3/cogent3",
    # turns off the secondary side-bar
    # it's default value is ["page-toc", "edit-this-page", "sourcelink"]
    "logo": {
        "text": "cogent3",
        "alt_text": "cogent3",
    },
    "collapse_navigation": False,
}

nbsphinx_thumbnails = make_nbsphinx_thumbnails()

# -- Options for LaTeX output --------------------------------------------------
latex_documents = [
    ("index", "cogent3.tex", "cogent3 Documentation", "cogent3 Team", "manual")
]
