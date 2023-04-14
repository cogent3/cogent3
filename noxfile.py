import pathlib

import nox


_py_versions = range(8, 11)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    session.install(".[test]")
    session.chdir("tests")
    session.run(
        "pytest",
        "-s",
        "-x",
        "--cov-report",
        f"lcov:lcov-{session.python}.info",
        "--cov",
        "cogent3",
        "--ignore",
        "test_app_mpi.py",
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
    for docdir in ("app", "cookbook", "examples"):
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
