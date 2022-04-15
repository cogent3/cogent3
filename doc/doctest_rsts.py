#!/usr/bin/env python
"""
This will doctest all files ending with .rst in this directory.
"""
import doctest
import os
import pathlib

import click
import nbformat

from nbconvert.preprocessors import CellExecutionError, ExecutePreprocessor

from cogent3.util.io import atomic_write

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


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


def execute_rsts(file_paths, exit_on_first, verbose):
    for test in file_paths:
        print()
        print("=" * 40)
        print(test)
        test = str(test)
        num_fails, num_tests = doctest.testfile(
            test,
            optionflags=doctest.ELLIPSIS or doctest.SKIP,
            verbose=verbose,
            encoding="utf-8",
        )
        if num_fails > 0 and exit_on_first:
            raise SystemExit(f"doctest failed in {test}")


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
        file_paths = [file_paths]
    else:
        file_paths = list(pathlib.Path(file_paths).glob(f"*.{suffix}"))

    if verbose:
        print(file_paths)

    if just:
        just = just.split(",")
        new = []
        for fn in file_paths:
            for sub_word in just:
                if sub_word in fn:
                    new.append(fn)
        file_paths = new
    elif exclude:
        exclude = exclude.split(",")
        new = []
        for fn in file_paths:
            keep = True
            for sub_word in exclude:
                if sub_word in fn:
                    keep = False
                    break
            if keep:
                new.append(fn)
        file_paths = new

    if verbose:
        print("File paths, after filtering: %s" % str(file_paths))

    if suffix == "rst":
        execute_rsts(file_paths, exit_on_first, verbose)
    else:
        execute_ipynb(file_paths, exit_on_first, verbose)


if __name__ == "__main__":
    main()
