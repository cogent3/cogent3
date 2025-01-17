import pathlib

import pytest


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
