import os
import shutil
import sys

from glob import glob

from sphinx_gallery.sorting import ExplicitOrder, FileNameSortKey

# set the plotly renderer
os.environ["PLOTLY_RENDERER"] = "sphinx_gallery"

sys.path.append("../src")

# Allow autosummary to generate stub files
autosummary_generate = True

# Prevent numpydoc from requiring stub files for methods
numpydoc_class_members_toctree = False

extensions = [
    "numpydoc",
    "sphinx.ext.todo",
    "sphinx.ext.doctest",
    "sphinx.ext.githubpages",
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
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "cookbook/union_dict.rst",
    "draw_examples/README.rst",
    "draw_examples/aln/README.rst",
    "draw_examples/tree/README.rst",
]

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
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
        # 'donate.html',
    ]
}


html_theme_options = {
    "fixed_sidebar": True,
}

html_static_path = ["_static"]

htmlhelp_basename = "cogent3doc"

# -- Options for Sphinx Gallery

nbsphinx_requirejs_path = "require.js"


def plotly_sg_scraper(block, block_vars, gallery_conf, *args, **kwargs):
    examples_dirs = gallery_conf["examples_dirs"]
    if isinstance(examples_dirs, (list, tuple)):
        examples_dirs = examples_dirs[0]
    pngs = sorted(glob(os.path.join(examples_dirs, "**", "*.png"), recursive=True))
    htmls = sorted(glob(os.path.join(examples_dirs, "**", "*.html"), recursive=True))
    image_path_iterator = block_vars["image_path_iterator"]
    image_names = list()
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
    try:
        result = plotly_sg_scraper(*args, **kwargs)
    except Exception:
        result = ""
    return result


image_scrapers = (plotly_sg_scraper_nb,)

example_dirs = ["draw_examples"]
gallery_dirs = ["draw"]


sphinx_gallery_conf = {
    # "doc_module": ("plotly",),
    "examples_dirs": example_dirs,
    "subsection_order": ExplicitOrder(["draw_examples/aln", "draw_examples/tree"]),
    "abort_on_example_error": True,
    "within_subsection_order": FileNameSortKey,
    "gallery_dirs": gallery_dirs,
    # "reference_url": {"plotly": None,
    # },
    "image_scrapers": image_scrapers,
    "download_all_examples": False,
}

# -- Options for LaTeX output --------------------------------------------------
latex_documents = [
    ("index", "cogent3.tex", "cogent3 Documentation", "cogent3 Team", "manual")
]
