#
# The following code is used for dynamic app installation and uninstallation in tests
#
# The general structure of an app to be dynamically loaded is:
#    @define_app
#    class MyAppClass:
#        def main(self, data: c3types.SerialisableType) -> str:
#            return json.dumps(data)
#        @classmethod
#        def _requires_imports(cls):  # imports required by the app
#            return ["from cogent3.app import typing as c3types", "import json"]
#
# Note: Due to the fact that it uses import_lib caches stevedore will see the first
# installed app immediately, but subsequent apps may not be seen until the cache is
# updated or the interpreter reloaded.


import importlib
import json
import os
import site
import sys
import time

from contextlib import contextmanager, suppress
from importlib.metadata import distributions
from inspect import getsourcelines, isclass
from os import path
from shutil import rmtree
from subprocess import check_call
from tempfile import mkdtemp

import pytest

import cogent3

from cogent3.app import apps, available_apps, get_app
from cogent3.app.composable import define_app


def is_package_installed(package_name):
    """Check if a package is installed"""
    site_packages_dir = site.getsitepackages()[1]
    installed_packages = os.listdir(site_packages_dir)
    return package_name in installed_packages


def install_app(cls, temp_dir, mod=None):
    """
    Installs a temporary app created from a class, by writing the source code
    into a new Python file in the temporary directory, creating a setup.py file
    that references the new module, and running pip install . the temporary directory.
    """
    # Get the source code of the function
    lines, _ = getsourcelines(cls)
    source_code = "\n".join(lines)
    if not mod:
        mod = f"mod_{cls.__name__}"

    # Write the source code into a new Python file in the temporary directory
    with open(path.join(temp_dir, f"{mod}.py"), "w") as f:
        f.write("from cogent3.app.composable import define_app\n")
        if hasattr(cls, "_requires_imports"):
            for import_statement in cls._requires_imports():
                f.write(f"{import_statement}\n")
        f.write(source_code)

    # Create a setup.py file in the temporary directory that references the new module
    setup_py = f"""
from setuptools import setup

setup(
name='{mod}',
version='0.1',
py_modules=['{mod}'],
entry_points={{
    "cogent3.app": [
        "{cls.__name__} = {mod}:{cls.__name__}",
    ],
}},    
)
"""
    with open(path.join(temp_dir, "setup.py"), "w") as f:
        f.write(setup_py)

    # Use subprocess to run pip install . in the temporary directory
    check_call([sys.executable, "-m", "pip", "install", "."], cwd=temp_dir)

    # Wait for the installation to complete
    timeout = 60  # maximum time to wait in seconds
    start_time = time.time()

    start_time = time.time()
    package_name = mod

    while True:
        if is_package_installed(package_name=package_name):
            break
        elif time.time() - start_time > timeout:
            raise TimeoutError(
                f"Package {package_name!r} not found after waiting for {timeout} seconds."
            )
        time.sleep(1)

    # check we can import it
    timeout = 60  # maximum time to wait in seconds
    start_time = time.time()

    while True:
        with suppress(ImportError):
            importlib.import_module(package_name)
            break
        if time.time() - start_time > timeout:
            raise TimeoutError(
                f"Package {package_name!r} not importable after waiting for {timeout} seconds."
            )
        time.sleep(1)

    # check it's in the stevedore cache
    timeout = 60  # maximum time to wait in seconds
    start_time = time.time()

    while True:
        appcache = apps(force=True)
        if any(ext.entry_point.value == f"{mod}:{cls.__name__}" for ext in appcache):
            break
        elif time.time() - start_time > timeout:
            raise TimeoutError(
                f"App {mod}.{cls.__name__} not available after waiting for {timeout} seconds."
            )
        time.sleep(1)


def uninstall_app(cls, mod=None):
    """
    Uninstalls a temporary app created from a class by running pip uninstall -y in the temporary directory.
    """
    if not mod:
        mod = f"mod_{cls.__name__}"
    # Use subprocess to run pip uninstall -y in the temporary directory
    package_name = mod.replace("_", "-")  # replace underscores with hyphens
    check_call([sys.executable, "-m", "pip", "uninstall", "-y", package_name])
    # force the app cache to reload
    _ = available_apps(force=True)


@contextmanager
def temp_app(cls, module_name: str = None):
    """
    A context manager that creates a temporary app from a class, installs it,
    and then uninstalls it after testing.

    Parameters
    ----------
    cls : object
        The class from which to create the app.
    module_name : str, optional
        The name of the module in which the app is defined. If None, a name
        is generated based on the class name.

    Yields
    ------
    None

    Raises
    ------
    pytest.fail
        If a TimeoutError occurs during app installation.

    Notes
    -----
    The app is installed in a temporary directory, which is deleted after testing.

    Usage
    -----
    with temp_app(cls):
        myapp = get_app(cls.__name__)
        # test myapp
    """
    if module_name is None:
        module_name = f"mod_{cls.__name__}"
    temp_dir = mkdtemp()
    try:
        try:
            install_app(cls, temp_dir=temp_dir, mod=module_name)
        except TimeoutError as e:
            pytest.fail(e.args[0])
        except Exception as e:
            pytest.fail(f"App installation failed: {e}")
        yield
    finally:
        uninstall_app(cls, mod=module_name)
        rmtree(temp_dir)


@pytest.fixture(scope="function")
def install_temp_app(request: pytest.FixtureRequest):
    """
    A pytest fixture that creates a temporary app from a class, installs it,
    and then uninstalls it after testing.

    Parameters
    ----------
    request : pytest.FixtureRequest
        The fixture request object. The class from which to create the app
        should be passed as a parameter to the test function using this fixture.

    Yields
    ------
    None

    Raises
    ------
    TypeError
        If the parameter passed to the test function is not a class.

    Notes
    -----
    The app is installed in a temporary directory, which is deleted after testing.

    Usage
    -----

    @pytest.mark.parametrize("install_temp_app", [MyAppClass], indirect=True)
    def test_function(install_temp_app):
        myapp = get_app('MyAppClass')
        # test myapp
    """
    cls = request.param
    if not isclass(cls):
        raise TypeError(f"Expected a class, but got {type(cls).__name__}")

    with temp_app(cls):
        yield


@pytest.mark.xfail(
    reason="This test is expected to fail due to a bug with changing from using pkg_resources to importlib"
)
def test_install_temp_app():
    @define_app
    class MyAppClass:
        def main(self, data: cogent3.app.typing.SerialisableType) -> str:
            return json.dumps(data)

        @classmethod
        def _requires_imports(cls):
            return ["from cogent3.app import typing as c3types", "import json"]

    with temp_app(MyAppClass):
        myapp = get_app("MyAppClass")
        data = {"key": "value", "number": 42, "list": [1, 2, 3], "boolean": True}
        assert (
            myapp(data)
            == '{"key": "value", "number": 42, "list": [1, 2, 3], "boolean": true}'
        )
