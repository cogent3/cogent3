import pathlib

from cogent3.util.io import (
    _path_relative_to_zip_parent,
)


def test_path_relative_to_zip_parent():
    """correctly generates member paths for a zip archive"""
    zip_path = pathlib.Path("some/path/to/a/data.zip")
    for member in ("data/member.txt", "member.txt", "a/b/c/member.txt"):
        got = _path_relative_to_zip_parent(zip_path, pathlib.Path(member))
        assert got.parts[0] == "data"
