import pathlib

import nox


dependencies = (
    "numba>0.54",
    "chardet",
    "numpy",
    "tinydb",
    "tqdm",
    "click",
    "pytest",
    "scitrack",
    "pandas",
    "plotly",
    "pytest-cov",
)

_py_versions = range(7, 11)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    py_version = session.python.replace(".", "")
    session.install(*dependencies)
    session.install(".")
    session.chdir("tests")
    session.run(
        "pytest",
        "-x",
        "--junitxml",
        f"junit-{py_version}.xml",
        "--cov-report",
        "xml",
        "--cov",
        "cogent3",
        "--ignore",
        "test_app_mpi.py",
    )


@nox.session(python=[f"3.{v}" for v in _py_versions])
def testmpi(session):
    session.install(*dependencies + ("mpi4py",))
    py = pathlib.Path(session.bin_paths[0]) / "python"
    session.install(".")
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
