import os

extensions = ['sphinx.ext.todo',
              'sphinx.ext.doctest',
#              'sphinx.ext.imgmath',
#              'nbsphinx',
              'sphinx.ext.mathjax',
#              'sphinxcontrib.bibtex'
              ]

# todo_include_todos=True # to expose the TODOs, uncomment this line

# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates']

# The suffix of source filenames.
source_suffix = '.rst'

# ignore the cookbook/ensembl.rst file as it's specifically imported
exclude_patterns = ['_build', '**.ipynb_checkpoints']

# The encoding of source files.
# source_encoding = 'utf-8'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'Cogent3'
copyright = u'2019, Cogent3'

version = ""

release = "2019.9.13a"

exclude_trees = ['_build']

show_authors = True

pygments_style = 'sphinx'

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  Major themes that come with
# Sphinx are currently 'default' and 'sphinxdoc'.

# on_rtd is whether we are on readthedocs.org, this line of code grabbed from docs.readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_static_path = ['_static']

htmlhelp_basename = 'Cogent3doc'


# -- Options for LaTeX output --------------------------------------------------
latex_documents = [
  ('index', 'Cogent3.tex', u'Cogent3 Documentation',
   u'Cogent3 Team', 'manual'),
]
