"""sets as working directory the parent of data/

Assumes this file is distributed inside COGENT3_ROOT/doc directory"""

import os
import pathlib

current = pathlib.Path(__file__).absolute().parent


def _data_path(path: os.PathLike) -> os.PathLike:
    if path.name != "data":
        path = _data_path(path.parent)
    return path


def get_data_dir():
    """returns path to cogent3 doc data directory"""
    for path in current.glob("*/*"):
        if "doc" not in path.parts:
            continue

        if "data" in path.parts:
            if path.parts.index("doc") > path.parts.index("doc"):
                continue

            path = _data_path(path)

        if path.is_dir() and str(path.name) == "data":
            return path.parent

    raise RuntimeError(f"could not find data dir from {current}")


def get_thumbnail_dir():
    """returns path to directory for writing html thumbnail images"""
    thumbdir = pathlib.Path(__file__).parent / "_build" / "html" / "_images"
    thumbdir.mkdir(exist_ok=True, parents=True)
    return thumbdir


data_dir = get_data_dir()
os.chdir(data_dir)
