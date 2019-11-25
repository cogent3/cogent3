#!/usr/bin/env python
import os
import re
import subprocess
import sys
import pathlib

from setuptools import Command, find_packages, setup
from setuptools.extension import Extension


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__contributors__ = [
    "Peter Maxwell",
    "Gavin Huttley",
    "Matthew Wakefield",
    "Greg Caporaso",
    "Daniel McDonald",
]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# Check Python version, no point installing if unsupported version inplace
min_version = (3, 6)
if sys.version_info < min_version:
    py_version = ".".join([str(n) for n in sys.version_info])
    msg = (
        f"Python-{'.'.join(min_version)} or greater is required, "
        f"Python-{py_version} used."
    )
    raise RuntimeError(msg)


# Check Numpy version, no point installing if unsupported version inplace
try:
    import numpy
except ImportError:
    raise RuntimeError("Numpy required but not found.")

numpy_version = re.split("[^\d]", numpy.__version__)
numpy_version_info = tuple([int(i) for i in numpy_version if i.isdigit()])
if numpy_version_info < (1, 3):
    raise RuntimeError("Numpy-1.3 is required, %s found." % numpy_version)

# Find arrayobject.h on any system
numpy_include_dir = numpy.get_include()


# On windows with no commandline probably means we want to build an installer.
if sys.platform == "win32" and len(sys.argv) < 2:
    sys.argv[1:] = ["bdist_wininst"]


# A new command for predist, ie: pyrexc but no compile.
class NullCommand(Command):
    description = "Generate .c files from .pyx files"
    # List of option tuples: long name, short name (or None), and help string.
    user_options = []  # [('', '', ""),]

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        pass


class BuildDocumentation(NullCommand):
    description = "Generate HTML documentation and .c files"

    def run(self):
        # Restructured Text -> HTML
        try:
            import sphinx
        except ImportError:
            print("Failed to build html due to ImportErrors for sphinx")
            return
        cwd = os.getcwd()
        os.chdir("doc")
        subprocess.call(["make", "html"])
        os.chdir(cwd)
        print("Built index.html")


# Cython is now run via the Cythonize function rather than monkeypatched into
# distutils, so these legacy commands don't need to do anything extra.
extra_commands = {
    "pyrexc": NullCommand,
    "cython": NullCommand,
    "predist": BuildDocumentation,
}


# Compiling Pyrex modules to .c and .so, if possible and necessary
try:
    if "DONT_USE_CYTHON" in os.environ:
        raise ImportError
    from Cython.Compiler.Version import version

    version = tuple([int(v) for v in re.split("[^\d]", version) if v.isdigit()])
    if version < (0, 17, 1):
        print("Your Cython version is too old")
        raise ImportError
except ImportError:
    source_suffix = ".c"
    cythonize = lambda x: x
    print("No Cython, will compile from .c files")
    for cmd in extra_commands:
        if cmd in sys.argv:
            print("'%s' command not available without Cython" % cmd)
            sys.exit(1)
else:
    from Cython.Build import cythonize

    source_suffix = ".pyx"


# Save some repetitive typing.  We have all compiled modules in place
# with their python siblings, so their paths and names are the same.
def CythonExtension(module_name, **kw):
    path = [PACKAGE_DIR] + module_name.split(".")
    path = os.path.join(*path) + source_suffix
    return Extension(module_name, [path], **kw)


short_description = "COmparative GENomics Toolkit 3"

readme_path = pathlib.Path(__file__).parent / "README.md"

long_description = readme_path.read_text()

PACKAGE_DIR = "src"

setup(
    name="cogent3",
    version=__version__,
    url="https://github.com/cogent3/cogent3",
    author="Gavin Huttley",
    author_email="gavin.huttley@anu.edu.au",
    description=short_description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    platforms=["any"],
    license=["BSD"],
    keywords=[
        "biology",
        "genomics",
        "statistics",
        "phylogeny",
        "evolution",
        "bioinformatics",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    packages=find_packages(where="src"),
    package_dir={"": PACKAGE_DIR},
    install_requires=["numpy", "pandas", "plotly", "scitrack", "tqdm", "tinydb"],
    extras_require={
        "dev": [
            "pytest-azurepipelines",
            "jupyterlab",
            "ipykernel",
            "ipywidgets",
            "click",
            "sphinx",
            "sphinx_rtd_theme",
            "sphinx-autobuild",
            "sphinxcontrib-bibtex",
            "numpydoc",
            "nbformat",
            "nbconvert",
            "jupyter_client",
            "ipykernel",
            "nbsphinx",
            "jupytext",
            "pytest",
            "pytest-cov",
            "pytest>=4.3.0",
            "tox",
            "black",
            "isort",
        ]
    },
    ext_modules=cythonize(
        [
            CythonExtension("cogent3.align._compare"),
            CythonExtension("cogent3.align._pairwise_seqs"),
            CythonExtension("cogent3.align._pairwise_pogs"),
            CythonExtension("cogent3.evolve._solved_models"),
            CythonExtension("cogent3.evolve._likelihood_tree"),
            CythonExtension("cogent3.evolve._pairwise_distance"),
            CythonExtension("cogent3.maths._period"),
        ]
    ),
    include_dirs=[numpy_include_dir],
    cmdclass=extra_commands,
)
