import pathlib

import nox

_py_versions = range(10, 13)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test_slow(session):
    session.install(".[test]")
    session.chdir("tests")
    session.run(
        "pytest",
        "-m",
        "slow",
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    session.install(".[test]")
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
        "--ignore",
        "test_app_mpi.py",
        "-m",
        "not slow",
        *session.posargs,
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test_module_docs(session):
    """doctest examples in a module"""
    session.install(".[test]")
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
    session.install(".[test]")
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
    session.install(".[doc]")
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


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test_new_type(session):
    session.install(".[test]")
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
