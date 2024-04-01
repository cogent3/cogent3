import importlib
import sys
import time

from contextlib import contextmanager
from importlib.metadata import EntryPoint, distribution
from inspect import getsourcelines, isclass
from os import path
from shutil import rmtree
from subprocess import check_call
from tempfile import mkdtemp
from unittest.mock import patch

import pytest

from stevedore import extension
from stevedore.extension import ExtensionManager

import cogent3

from cogent3.app import (
    _make_apphelp_docstring,
    app_help,
    apps,
    available_apps,
    get_app,
)
from cogent3.app.composable import define_app
from cogent3.util.table import Table


@pytest.fixture
def extension_manager_factory():
    """Fixture to create mocked ExtensionManager instances with given extensions."""

    def _factory(extensions):
        with patch("stevedore.ExtensionManager") as mock_manager_constructor:
            mock_manager = ExtensionManager.make_test_instance(
                extensions=extensions, namespace="TESTING"
            )
            mock_manager_constructor.return_value = mock_manager
            return mock_manager

    return _factory


@pytest.fixture
def mock_extension_manager(extension_manager_factory, monkeypatch):
    """Fixture to mock the ExtensionManager with the given extensions."""

    def _mock_extension_manager(extensions):
        # Create a mocked ExtensionManager with the specified mock extensions
        mocked_manager = extension_manager_factory(extensions)
        # Patch the __apps variable in cogent3.apps module to use the mocked_manager
        monkeypatch.setattr(cogent3.app, "__apps", mocked_manager)
        return mocked_manager

    return _mock_extension_manager


def create_extension(
    plugin: object, name: str = None, module_name: str = "module1"
) -> extension.Extension:
    if name is None:
        name = plugin.__name__
    return extension.Extension(
        name=name,
        entry_point=EntryPoint(
            name=name, value=f"{module_name}:{name}", group="TESTING"
        ),
        plugin=plugin,
        obj=None,
    )


def test_Install_app_class(mock_extension_manager):
    @define_app
    class uppercase:
        """Test app that converts a string to uppercase"""

        def main(self, data: str) -> str:
            return data.upper()

    mock_extension_manager([create_extension(uppercase)])

    appercase = get_app("uppercase")
    assert appercase("hello") == "HELLO"
    assert appercase.__doc__ in _make_apphelp_docstring(appercase.__class__)


def test_Install_app_function(mock_extension_manager):
    @define_app
    def uppercase(data: str) -> str:
        """Test function that converts a string to uppercase"""
        return data.upper()

    mock_extension_manager([create_extension(uppercase)])

    appercase = get_app("uppercase")
    assert appercase("hello") == "HELLO"
    assert appercase.__doc__ in _make_apphelp_docstring(appercase.__class__)


@pytest.mark.parametrize("app_doc", [None, "text"])
@pytest.mark.parametrize("init_doc", [None, "text"])
def test_app_docs(mock_extension_manager, app_doc, init_doc):
    @define_app
    class documented_app:
        """This is a test app that has a __init__, and a docstring"""

        def __init__(self):
            self.constant = 2

        def main(self, val: int) -> int:
            return val + self.constant

    mock_extension_manager([create_extension(documented_app)])

    assert cogent3.app.apps().names() == ["documented_app"]
    app = get_app("documented_app")
    app.__class__.__doc__ = app_doc
    app.__class__.__init__.__doc__ = init_doc
    app_help("documented_app")
    got = _make_apphelp_docstring(app.__class__)
    assert "Options" in got


def test_namespace_collision(mock_extension_manager):
    @define_app
    class app1:
        def main(self, data: str) -> str:
            return data.upper()

    @define_app
    class app2:
        def main(self, data: str) -> str:
            return data.lower()

    # create two apps with the same name and different modules
    mock_extension_manager(
        [
            create_extension(app1, module_name="module1"),
            create_extension(app2, name="app1", module_name="module2"),
        ]
    )

    assert cogent3.app.apps().names() == ["app1", "app1"]

    with pytest.raises(NameError):
        _ = get_app(
            "app1"
        )  # request app by name only, when there are multiple apps with the same name

    app_by_module_name_1 = get_app("module1.app1")  # request app by name and module
    app_by_module_name_2 = get_app("module2.app1")
    assert app_by_module_name_1("Hello") == "HELLO"
    assert app_by_module_name_2("Hello") == "hello"

    composition = app_by_module_name_1 + app_by_module_name_2
    assert composition("Hello") == "hello"


def test_available_apps_local(mock_extension_manager):
    """available_apps robust to local scope apps"""

    @define_app
    def dummy(val: int) -> int:
        return val

    mock_extension_manager([create_extension(dummy)])
    apps = available_apps()
    assert isinstance(apps, Table)
    apps.filtered(lambda x: dummy.__name__ == x, columns="name")
    assert apps.shape[0] == 1


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
    while True:
        importlib.invalidate_caches()
        installed_packages = [dist.metadata["Name"] for dist in distribution()]
        package_name = mod.replace("_", "-")  # replace underscores with hyphens
        if package_name in installed_packages:
            print(f"Package {package_name!r} found.")
            break
        elif time.time() - start_time > timeout:
            raise TimeoutError(
                f"Package {package_name!r} not found after waiting for {timeout} seconds."
            )
        time.sleep(1)

    # wait until the stevedore cache is updated
    timeout = 60  # maximum time to wait in seconds
    start_time = time.time()

    while True:
        appcache = apps(force=True)
        if any(ext.entry_point.value == f"{mod}:{cls.__name__}" for ext in appcache):
            break
        elif time.time() - start_time > timeout:
            raise TimeoutError(
                f"App {mod}.{cls.__name__} not found after waiting for {timeout} seconds."
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
    apps = available_apps(force=True)


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
        except TimeoutError:
            pytest.fail(
                f"TimeoutError occurred during `{cls.__name__}` app installation"
            )
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
