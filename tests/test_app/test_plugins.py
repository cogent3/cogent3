import importlib
import json
import sys
import time

from inspect import getsourcelines, isclass
from os import path
from shutil import rmtree
from subprocess import check_call, check_output
from tempfile import mkdtemp

import pkg_resources
import pytest

from cogent3.app import _make_apphelp_docstring, app_help, available_apps, get_app
from cogent3.app import typing as c3types
from cogent3.app.composable import define_app


@pytest.fixture(scope="function")
def install_temp_app(request: pytest.FixtureRequest):
    """
    A pytest fixture that creates a temporary app from a function, installs it, and then uninstalls it after testing.

    This fixture creates a temporary directory, writes the source code of the function into a new Python file in the
    temporary directory, creates a setup.py file that references the new module, and installs the module using pip.
    After the test function that uses this fixture is run, the module is uninstalled and the temporary directory is deleted.

    Parameters
    ----------
    request : _pytest.fixtures.SubRequest
        An define_app decorated class

    Notes
    -----

    This fixture will only work with define_app decorated classes, NOT define_app decorated functions that have been shoehorned
    into classes because with a dynamically generated class the original code is not available for inspection programatically.

    If your temporary app will require special imports, you can add import statements in a class method called `_requires_imports`.
    This method should return a list of strings that are valid import statements.

    Yields
    ------
    None

    Raises
    ------
    subprocess.CalledProcessError
        If the `pip install .` or `pip uninstall -y` commands fail.

    Side Effects
    ------------
    - forces refresh of the apps cache

    Examples
    --------
    >>> @pytest.mark.parametrize('install_temp_module', [my_func], indirect=True)
    ... def test_my_func(install_temp_module):
    ...     # Test code here
    """
    # Create a temporary directory
    temp_dir = mkdtemp()

    # Get the function from the request
    cls = request.param

    # Check if the passed-in value is a class ( a define_app decorated function becomes a class)
    if not isclass(cls):
        raise TypeError(f"Expected a class, but got {type(cls).__name__}")

    # Get the source code of the function
    lines, _ = getsourcelines(cls)
    source_code = "\n".join(lines)

    # Write the source code into a new Python file in the temporary directory
    with open(path.join(temp_dir, f"{cls.__name__}.py"), "w") as f:
        f.write("from cogent3.app.composable import define_app\n")
        if hasattr(cls, "_requires_imports"):
            for import_statement in cls._requires_imports():
                f.write(f"{import_statement}\n")
        f.write(source_code)

    # Create a setup.py file in the temporary directory that references the new module
    setup_py = f"""
from setuptools import setup

setup(
    name='{cls.__name__}',
    version='0.1',
    py_modules=['{cls.__name__}'],
    entry_points={{
        "cogent3.app": [
            "{cls.__name__} = {cls.__name__}:{cls.__name__}",
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

    while True:
        installed_packages = check_output(
            [sys.executable, "-m", "pip", "list"]
        ).decode()
        package_name = cls.__name__.replace(
            "_", "-"
        )  # replace underscores with hyphens
        if package_name in installed_packages:
            break
        elif time.time() - start_time > timeout:
            raise TimeoutError(
                f"Package {package_name} not found after waiting for {timeout} seconds."
            )
        time.sleep(1)

    # force pkg_resourcs to reload it's cache of entry points
    importlib.reload(pkg_resources)

    # Create a new ExtensionManager instance to force it to reload the entry points
    em = ExtensionManager(namespace="cogent3.app", invoke_on_load=False)

    # force the app cache to reload
    apps = available_apps(force=True)

    yield

    # After the test function completes, use subprocess to run pip uninstall -y <module> to uninstall the module
    check_call([sys.executable, "-m", "pip", "uninstall", "-y", f"{cls.__name__}"])

    # Create a new ExtensionManager instance to force it to reload the entry points
    em = ExtensionManager(namespace="cogent3.app", invoke_on_load=False)

    # force the app cache to reload after
    apps = available_apps(force=True)

    # Delete the temporary directory
    rmtree(temp_dir)


@define_app
class uppercase:
    """Convert a string to uppercase."""
    def main(self, data: str) -> str:
        return data.upper()


@pytest.mark.parametrize("install_temp_app", [uppercase], indirect=True)
def test_add_app(install_temp_app):
    appercase = get_app("uppercase")
    assert appercase("hello") == "HELLO"
    assert appercase.__doc__ in _make_apphelp_docstring(appercase.__class__)


@define_app
class to_json:
    def main(self, data: c3types.SerialisableType) -> str:
        """Convert primitive python types to json string."""
        return json.dumps(data)

    @classmethod
    def _requires_imports(cls):
        return ["from cogent3.app import typing as c3types"]


from stevedore.extension import ExtensionManager


@pytest.mark.parametrize("install_temp_app", [to_json], indirect=True)
def test_app_namespace_collision(install_temp_app):
    pr = pkg_resources.iter_entry_points(group="cogent3.app")
    to_json_extensions = [
        entry_point for entry_point in pr if entry_point.name == "to_json"
    ]
    assert (
        len(to_json_extensions) == 2
    ), "pkg_resources.entry_points does not show 2 to_json apps"

    em = ExtensionManager(namespace="cogent3.app", invoke_on_load=False)
    to_json_extensions = [ext for ext in em.extensions if ext.name == "to_json"]
    assert (
        len(to_json_extensions) == 2
    ), "stevedore.ExtensionManager does not show 2 to_json apps"

    with pytest.raises(NameError):
        to_json = get_app("to_json")

    assert (
        available_apps("to_json").shape[0] == 2
    ), "available_apps does not show 2 to_json apps"
