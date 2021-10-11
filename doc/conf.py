import datetime
import os
import shutil
import sys

from glob import glob

import sphinx_bootstrap_theme

from sphinx_gallery.sorting import ExplicitOrder, FileNameSortKey


sys.path.append("../src")

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
    "sphinx_gallery.gen_gallery",
    "sphinx_panels",
    "sphinxcontrib.bibtex",
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
    "draw_examples/README.rst",
    "draw_examples/aln/README.rst",
    "draw_examples/tree/README.rst",
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

release = "2021.10.12a"

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

# -- Options for Sphinx Gallery

nbsphinx_requirejs_path = "require.js"

# This is a veryt ugly hack!
# We have to set the following environment variable
# for Sphinx-Gallery to work, but we need the default Plotly setting for correct
# display in normal notebooks. This works because of an order of execution,
# over which we have no control, putting the gallery generation first.
os.environ["PLOTLY_RENDERER"] = "sphinx_gallery"


def plotly_sg_scraper(block, block_vars, gallery_conf, *args, **kwargs):
    examples_dirs = gallery_conf["examples_dirs"]
    if isinstance(examples_dirs, (list, tuple)):
        examples_dirs = examples_dirs[0]
    pngs = sorted(glob(os.path.join(examples_dirs, "**", "*.png"), recursive=True))
    htmls = sorted(glob(os.path.join(examples_dirs, "**", "*.html"), recursive=True))
    image_path_iterator = block_vars["image_path_iterator"]
    image_names = []
    seen = set()
    for html, png in zip(htmls, pngs):
        if png not in seen:
            seen |= {png}
            this_image_path_png = next(image_path_iterator)
            this_image_path_html = os.path.splitext(this_image_path_png)[0] + ".html"
            image_names.append(this_image_path_html)
            shutil.move(png, this_image_path_png)
            shutil.move(html, this_image_path_html)
    # Use the `figure_rst` helper function to generate rST for image files
    from plotly.io._sg_scraper import figure_rst

    return figure_rst(image_names, gallery_conf["src_dir"])


def plotly_sg_scraper_nb(*args, **kwargs):
    # set the plotly renderer
    os.environ["PLOTLY_RENDERER"] = "sphinx_gallery"
    try:
        result = plotly_sg_scraper(*args, **kwargs)
    except Exception:
        result = ""

    os.environ.pop("PLOTLY_RENDERER")
    return result


image_scrapers = (plotly_sg_scraper_nb,)

example_dirs = ["draw_examples"]
gallery_dirs = ["draw"]


sphinx_gallery_conf = {
    "examples_dirs": example_dirs,
    "subsection_order": ExplicitOrder(["draw_examples/aln", "draw_examples/tree"]),
    "abort_on_example_error": True,
    "within_subsection_order": FileNameSortKey,
    "gallery_dirs": gallery_dirs,
    "image_scrapers": image_scrapers,
    "download_all_examples": False,
}

# -- Options for LaTeX output --------------------------------------------------
latex_documents = [
    ("index", "cogent3.tex", "cogent3 Documentation", "cogent3 Team", "manual")
]


def setup(app):
    app.add_js_file("plotly-latest.min.js")
