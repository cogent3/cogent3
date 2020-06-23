"""sets as working directory the parent of data/

Assumes this file is distributed inside COGENT3_ROOT/doc directory"""
import os
import pathlib


def get_data_dir():
    """returns path to cogent3 doc data directory"""
    current = pathlib.Path(".").absolute().parent
    for path in current.glob("**/*"):
        if "doc" not in path.parts:
            continue

        if path.is_dir() and str(path.name) == "data":
            return path.parent

    raise RuntimeError(f"could not find data dir from {current}")


data_dir = get_data_dir()
os.chdir(data_dir)
