import os
import pathlib

import pytest


def pytest_configure() -> None:
    # tests use paths relative to this directory
    os.chdir(pathlib.Path(__file__).parent)


@pytest.fixture(scope="session")
def DATA_DIR() -> pathlib.Path:
    return pathlib.Path(__file__).parent / "data"


@pytest.fixture
def HOME_TMP_DIR(DATA_DIR) -> pathlib.Path:
    """makes a temporary directory"""
    import tempfile

    HOME = pathlib.Path("~")
    with tempfile.TemporaryDirectory(dir=HOME.expanduser()) as dn:
        yield HOME / dn
