import os
import pathlib
import sys

import nox

# on python >= 3.12 this will improve speed of test coverage a lot
if sys.version_info >= (3, 12):
    os.environ["COVERAGE_CORE"] = "sysmon"

_py_versions = range(11, 16)

_typed_modules = [
    "src/cogent3/core/alignment.py",
    "src/cogent3/core/alphabet.py",
    "src/cogent3/core/annotation.py",
    "src/cogent3/core/genetic_code.py",
    "src/cogent3/core/location.py",
    "src/cogent3/core/moltype.py",
    "src/cogent3/core/seq_storage.py",
    "src/cogent3/core/sequence.py",
    "src/cogent3/core/seqview.py",
    "src/cogent3/core/slice_record.py",
    "src/cogent3/core/tree.py",
]


@nox.session(python=False)
def fmt(session: nox.Session) -> None:
    session.run("ruff", "check", "--fix-only", ".", external=True)
    session.run("ruff", "format", ".", external=True)


@nox.session(python=[f"3.{v}" for v in _py_versions], venv_backend="uv")
def test_slow(session):
    session.install("-e", ".", "--group", "test")
    session.chdir("tests")
    session.run(
        "pytest",
        "-m",
        "slow",
    )


@nox.session(python=[f"3.{v}" for v in _py_versions], venv_backend="uv")
def test(session):
    session.install("-e", ".", "--group", "test")
    session.run("uv", "pip", "list", external=True)
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


@nox.session(python="3.14", venv_backend="uv")
def type_check(session):
    session.install("-e", ".", "--group", "typing")
    session.run("mypy", *_typed_modules)


@nox.session(python=[f"3.{v}" for v in _py_versions], venv_backend="uv")
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


@nox.session(python=[f"3.{v}" for v in _py_versions], venv_backend="uv")
def testdocs(session):
    """render the docs, which executes every code cell"""
    py = pathlib.Path(session.bin_paths[0]) / "python"
    session.install("-e", ".", "--group", "doc")
    # quarto discovers the interpreter, and so the ipykernel to execute cells
    # with, from QUARTO_PYTHON
    session.env["QUARTO_PYTHON"] = str(py)
    session.run("quartodoc", "build", "--config", "doc/_quarto.yml")
    session.run("quarto", "render", "doc", *session.posargs, external=True)
