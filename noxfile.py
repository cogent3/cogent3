import os
import pathlib
import sys

import nox

# on python >= 3.12 this will improve speed of test coverage a lot
if sys.version_info >= (3, 12):
    os.environ["COVERAGE_CORE"] = "sysmon"

_py_versions = range(11, 15)


@nox.session(python=False)
def fmt(session: nox.Session) -> None:
    session.run("ruff", "check", "--fix-only", ".", external=True)
    session.run("ruff", "format", ".", external=True)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test_slow(session):
    session.install("-e", ".", "--group", "test")
    session.chdir("tests")
    session.run(
        "pytest",
        "-m",
        "slow",
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    session.install("-e", ".", "--group", "test")
    session.run("pip", "list")
    # doctest modules within cogent3/app
    session.chdir("src/cogent3/app")
    session.run(
        "pytest",
        "-s",
        "-x",
        "--doctest-modules",
        "--ignore=composable.py",
        "--ignore=sqlite_data_store.py",
        ".",
    )

    session.chdir("../../../tests")
    session.run(
        "pytest",
        "-s",
        "-x",
        "-m",
        "not slow",
        *session.posargs,
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test_module_docs(session):
    """doctest examples in a module"""
    session.install("-e", ".", "--group", "test")
    # doctest modules within cogent3/app
    session.chdir("src/cogent3/app")
    session.run(
        "pytest",
        "-s",
        "--doctest-modules",
        "--ignore=composable.py",
        "--ignore=sqlite_data_store.py",
        ".",
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def testdocs(session):
    """render the docs, which executes every code cell"""
    py = pathlib.Path(session.bin_paths[0]) / "python"
    session.install("-e", ".", "--group", "doc")
    # quarto discovers the interpreter, and so the ipykernel to execute cells
    # with, from QUARTO_PYTHON
    session.env["QUARTO_PYTHON"] = str(py)
    session.run("quartodoc", "build", "--config", "doc/_quarto.yml")
    session.run("quarto", "render", "doc", *session.posargs, external=True)
