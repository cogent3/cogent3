import pathlib

import nox


_py_versions = range(7, 11)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    py_version = session.python.replace(".", "")
    session.install(".[test]")
    session.chdir("tests")
    session.run(
        "pytest",
        "-s",
        "-x",
        "--junitxml",
        f"junit-{session.python}.xml",
        "--cov-report",
        "xml",
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
    session.run(
        str(py),
        "doctest_rsts.py",
        "-f",
        "cookbook",
        "-1",
        "-s",
        "rst",
        external=True,
    )
    session.run(
        str(py),
        "doctest_rsts.py",
        "-f",
        "examples",
        "-1",
        "-s",
        "rst",
        external=True,
    )
