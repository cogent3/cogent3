import gc
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


@pytest.fixture(scope="session", autouse=True)
def _try_cleaning_up_on_autouse_fixture_teardown():
    yield
    for _ in range(10):
        gc.collect()


def pytest_sessionfinish(session, exitstatus):
    for _ in range(10):
        gc.collect()
