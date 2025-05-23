import os
import pathlib

import nox

# on python >= 3.12 this will improve speed of test coverage a lot
os.environ["COVERAGE_CORE"] = "sysmon"

_py_versions = range(10, 14)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test_slow(session):
    session.install("-e.[test]")
    session.chdir("tests")
    session.run(
        "pytest",
        "-m",
        "slow",
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    session.install("-e.[test]")
    session.run("pip", "list")
    # doctest modules within cogent3/app
    session.chdir("src/cogent3/app")
    session.run(
        "pytest",
        "-s",
        "-x",
        "--doctest-modules",
        ".",
    )

    session.chdir("../../../tests")
    session.run(
        "pytest",
        "-s",
        "-x",
        "--ignore-glob",
        "*app/test_app_mpi.py",
        "--ignore-glob",
        "*core/test_alignment.py",
        "--ignore-glob",
        "*core/test_alphabet.py",
        "--ignore-glob",
        "*core/test_genetic_code.py",
        "--ignore-glob",
        "*core/test_moltype.py",
        "--ignore-glob",
        "*core/test_sequence.py",
        "-m",
        "not slow",
        *session.posargs,
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test_module_docs(session):
    """doctest examples in a module"""
    session.install("-e.[test]")
    # doctest modules within cogent3/app
    session.chdir("src/cogent3/app")
    session.run(
        "pytest",
        "-s",
        "--doctest-modules",
        ".",
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def testmpi(session):
    session.install("-e.[test]")
    session.install("mpi4py")
    py = pathlib.Path(session.bin_paths[0]) / "python"
    session.chdir("tests")
    session.run(
        "mpiexec",
        "-n",
        "2",
        str(py),
        "-m",
        "mpi4py.futures",
        "-m",
        "pytest",
        "-x",
        "test_app/test_app_mpi.py",
        external=True,
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def testdocs(session):
    py = pathlib.Path(session.bin_paths[0]) / "python"
    session.install("-e.[doc]")
    session.chdir("doc")
    for docdir in ("app", "cookbook", "draw", "examples"):
        session.run(
            str(py),
            "doctest_rsts.py",
            "-f",
            docdir,
            "-1",
            "-s",
            "rst",
            external=True,
        )


@nox.session(python=None)
def test_new_type(session):
    session.env["COGENT3_NEW_TYPE"] = "1"
    session.chdir("tests")
    session.run(
        "pytest",
        "-s",
        "--ignore",
        "test_app_mpi.py",
        "--ignore",
        "test_core/test_alignment.py",
        "--ignore",
        "test_core/test_alphabet.py",
        "--ignore",
        "test_core/test_annotation.py",
        "--ignore",
        "test_core/test_genetic_code.py",
        "--ignore",
        "test_core/test_moltype.py",
        "--ignore",
        "test_core/test_sequence.py",
        "--ignore",
        "test_core/test_seq_aln_integration.py",
        "--ignore",
        "test_core/test_core_standalone.py",
        "--ignore",
        "test_util/test_recode_alignment.py",
        "-m",
        "not slow",
        *session.posargs,
    )
