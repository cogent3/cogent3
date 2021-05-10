import pathlib
import sys

from setuptools import find_packages, setup


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2021, The Cogent Project"
__contributors__ = [
    "Peter Maxwell",
    "Gavin Huttley",
    "Matthew Wakefield",
    "Greg Caporaso",
    "Daniel McDonald",
]
__license__ = "BSD-3"
__version__ = "2021.5.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

# Check Python version, no point installing if unsupported version inplace
min_version = (3, 7)
if sys.version_info < min_version:
    py_version = ".".join(str(n) for n in sys.version_info)
    msg = (
        f"Python-{'.'.join(map(str, min_version))} or greater is required, "
        f"Python-{py_version} used."
    )
    raise RuntimeError(msg)


# On windows with no commandline probably means we want to build an installer.
if sys.platform == "win32" and len(sys.argv) < 2:
    sys.argv[1:] = ["bdist_wininst"]


short_description = "COmparative GENomics Toolkit 3"

readme_path = pathlib.Path(__file__).parent / "README.md"

long_description = readme_path.read_text()


PACKAGE_DIR = "src"

PROJECT_URLS = {
    "Documentation": "https://www.cogent3.org/",
    "Bug Tracker": "https://github.com/cogent3/cogent3/issues",
    "Source Code": "https://github.com/cogent3/cogent3",
}

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
    license=__license__,
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
        "License :: OSI Approved :: BSD License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    packages=find_packages(where="src"),
    package_dir={"": PACKAGE_DIR},
    install_requires=[
        "chardet",
        "numpy",
        "numba>0.48.0;python_version<'3.9'",
        "numba>0.53; python_version>='3.9'",
        "scitrack",
        "tqdm",
        "tinydb",
    ],
    extras_require={
        "dev": [
            "black",
            "click",
            "ipykernel",
            "ipywidgets",
            "isort",
            "jupyter-sphinx",
            "jupyter_client",
            "jupyterlab",
            "jupytext",
            "kaleido",
            "nbconvert",
            "nbformat",
            "nbsphinx",
            "numpydoc",
            "pandas",
            "plotly",
            "psutil",
            "pytest",
            "pytest-cov",
            "pytest>=4.3.0",
            "sphinx",
            "sphinx-autobuild",
            "sphinxcontrib-bibtex",
            "sphinx_panels",
            "tox",
        ],
        "extra": ["pandas", "plotly", "psutil", "kaleido"],
    },
    project_urls=PROJECT_URLS,
)
