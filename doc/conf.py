import os
import pathlib
import shutil
import subprocess
import sys

from glob import glob

import plotly.io as pio

from sphinx_gallery.sorting import ExplicitOrder, FileNameSortKey

# following functions used to hack readthedocs build environment to include
# required tools for plotly to work on a headless linux box as root
def exec_command(cmnd):
    """executes shell command and returns stdout if completes exit code 0"""
    proc = subprocess.Popen(
        cmnd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = proc.communicate()
    if proc.returncode != 0:
        user = os.environ.get("USER", None)
        err = err.decode("utf8")
        raise SystemError(proc.returncode, f"FAILED: {cmnd}\n{err} for '{user}'")

    if out is not None:
        r = out.decode("utf8")
    else:
        r = None

    return r


def not_installed_on_linux(packages):
    # returns packages to be installed, if not Linux just an empty list
    if "linux" not in sys.platform.lower():
        return []

    to_install = {}
    for name, url in packages.items():
        cmnd = f"dpkg -l | grep {name}"
        try:
            result = exec_command(cmnd)
        except SystemError:
            to_install[name] = url

    return to_install


def update_orca():
    # write out a revised plot-orca command
    orca = pathlib.Path(exec_command("which orca").strip())
    txt = orca.read_text()
    # get the linux line
    cmnd = [l.strip() for l in txt.split("linux")[-1].splitlines() if "exec" in l]
    cmnd = f"{cmnd[0]} --no-sandbox"
    cmnd = [
        "#!/bin/bash",
        "",
        cmnd,
    ]
    # copy the original
    shutil.copy(str(orca), f"{orca}.orig")
    # write the new
    orca.write_text("\n".join(cmnd))


def apt_get_installs():
    # need this to get around issues of no X11 on readthedocs
    # not really doing an apt-get install, but an apt-get download of a .deb
    # then using dpkg to expand the packages locally
    if "linux" not in sys.platform.lower():
        return

    packages = not_installed_on_linux(
        {
            "libgtk2.0-0": "http://nova.clouds.archive.ubuntu.com/ubuntu/pool/main/g/gtk+2.0/libgtk2.0-0_2.24.30-1ubuntu1.16.04.2_amd64.deb",
            "libgconf-2-4": "http://nova.clouds.archive.ubuntu.com/ubuntu/pool/main/g/gconf/libgconf-2-4_3.2.6-3ubuntu6_amd64.deb",
            "xvfb": "http://nova.clouds.archive.ubuntu.com/ubuntu/pool/universe/x/xorg-server/xvfb_1.18.4-0ubuntu0.8_amd64.deb",
            "chromium-browser": "http://nova.clouds.archive.ubuntu.com/ubuntu/pool/universe/c/chromium-browser/chromium-browser_80.0.3987.87-0ubuntu0.16.04.1_amd64.deb",
        }
    )
    if not packages:
        return

    for name, package in packages.items():
        print(f"Installing {name}")
        cmnd = f"wget {package}"
        exec_command(cmnd)
        package = package.split("/")[-1]
        cmnd = f"dpkg -x {package} ."
        exec_command(cmnd)

    if packages:
        # now delete the *.deb
        exec_command("rm *.deb")

    # now modify PATH and LD_LIBRARY_PATH
    install_dir = os.getcwd()
    for env, local in [("PATH", "bin"), ("LD_LIBRARY_PATH", "lib")]:
        local = os.path.join(install_dir, "usr", local)
        e = os.environ.get(env, "")
        if local in e:
            continue

        os.environ[env] = ":".join([local, e])

    # modify orca for root exec on linux
    update_orca()


apt_get_installs()


try:
    # if Xvfb installed
    r = exec_command("which xvfb-run")
    pio.orca.config.use_xvfb = True
except SystemError as msg:
    if "linux" in sys.platform.lower():
        print(msg)

# end of hack attempt

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
