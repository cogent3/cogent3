#!/usr/bin/env python
"""
This will doctest all files ending with .rst in this directory.
"""
import os
import pathlib
import subprocess
import tempfile

import click
import nbformat

from nbconvert.preprocessors import CellExecutionError, ExecutePreprocessor

from cogent3.app.composable import define_app
from cogent3.util.io import atomic_write


DATA_DIR = pathlib.Path(__file__).parent / "data"


def execute_ipynb(file_paths, exit_on_first, verbose):
    failed = False
    for test in file_paths:
        path = pathlib.Path(test)
        print()
        print("=" * 40)
        print(test)
        with open(test) as f:
            nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(
            timeout=600, kernel_name="python3", store_widget_state=True
        )
        try:
            ep.preprocess(nb, {"metadata": {"path": path.parent}})
        except CellExecutionError:
            failed = True

        with atomic_write(test, mode="w") as f:
            nbformat.write(nb, f)

        if failed and exit_on_first:
            raise SystemExit(
                f"notebook execution failed in {test}, error saved " "in notebook"
            )


@define_app
def test_file(test: pathlib.Path | str, exit_on_first: bool = True) -> bool:
    test = pathlib.Path(test).absolute()
    orig_wd = pathlib.Path(os.getcwd())
    rst_to_py = orig_wd / "rst2script.py"
    with tempfile.TemporaryDirectory(dir=orig_wd) as working_dir:
        working_dir = pathlib.Path(working_dir)
        # create a symlink to the data dir
        (working_dir / "data").symlink_to(DATA_DIR)

        # copy the rest file
        dest = working_dir / test.name
        dest.write_text(test.read_text())
        test = dest
        cmnd = f"python {str(rst_to_py)} -wd {str(working_dir)} {str(test)}"
        r = subprocess.run(cmnd.split(), capture_output=True, check=False)
        if r.returncode != 0:
            click.secho(f"FAILED: {str(rst_to_py)}", fg="red")
            click.secho(r.stdout.decode("utf8"), fg="red")
            click.secho(r.stderr.decode("utf8"), fg="red")
            exit(1)

        py_path = pathlib.Path(f'{str(test).removesuffix("rst")}py')

        cmnd = f"python {str(py_path)}"
        r = subprocess.run(cmnd.split(), capture_output=True)

        if r.returncode != 0:
            click.secho(f"FAILED: {str(py_path)}", fg="red")
            click.secho(r.stdout.decode("utf8"), fg="red")
            click.secho(r.stderr.decode("utf8"), fg="red")

        if exit_on_first and r.returncode != 0:
            return False

    return True


def execute_rsts(file_paths, exit_on_first, verbose):
    runfile = test_file(exit_on_first=exit_on_first)
    for result in runfile.as_completed(file_paths, parallel=True):
        if result is False and exit_on_first:
            exit(1)


@click.command()
@click.option(
    "-f",
    "--file_paths",
    required=True,
    help="directory or specific"
    " files to test. If directory, glob searches for files matching suffix.",
)
@click.option(
    "-j",
    "--just",
    help="comma separated list of names to be matched to files to be tested",
)
@click.option(
    "-x",
    "--exclude",
    help="comma separated list of names to be matched to files to be excluded",
)
@click.option("-1", "--exit_on_first", is_flag=True, help="exit on first failure")
@click.option(
    "-s", "--suffix", type=click.Choice(["rst", "ipynb"]), help="suffix of docs to test"
)
@click.option("-v", "--verbose", is_flag=True, help="verbose output")
def main(file_paths, just, exclude, exit_on_first, suffix, verbose):
    """runs doctests for the indicated files"""
    cwd = os.getcwd()
    if "*" in file_paths:  # trim to just parent directory
        file_paths = pathlib.Path(file_paths).parent

    if "," in file_paths:
        file_paths = file_paths.split(",")
    elif pathlib.Path(file_paths).is_file():
        p = pathlib.Path(file_paths)
        assert p.suffix.lower() in ".rst.ipynb", f"Unknown format suffix '{p.suffix}'"
        suffix = p.suffix[1:]
        file_paths = [p]
    else:
        file_paths = list(pathlib.Path(file_paths).glob(f"*.{suffix}"))

    if verbose:
        print(file_paths)

    if just:
        just = just.split(",")
        new = []
        for fn in file_paths:
            new.extend(fn for sub_word in just if sub_word in fn.name)
        file_paths = new
    elif exclude:
        exclude = exclude.split(",")
        new = []
        for fn in file_paths:
            keep = all(sub_word not in str(fn) for sub_word in exclude)
            if keep:
                new.append(fn)
        file_paths = new

    if verbose:
        print(f"File paths, after filtering: {str(file_paths)}")

    if suffix == "rst":
        execute_rsts(file_paths, exit_on_first, verbose)
    else:
        execute_ipynb(file_paths, exit_on_first, verbose)


if __name__ == "__main__":
    main()
