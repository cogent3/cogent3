"""sets as working directory the parent of data/

Assumes this file is distributed inside COGENT3_ROOT/doc directory"""
import os
import pathlib


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def get_data_dir():
    """returns path to cogent3 doc data directory"""
    current = pathlib.Path(".").absolute().parent
    for path in current.glob("**/*"):
        if "doc" not in path.parts:
            continue

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
