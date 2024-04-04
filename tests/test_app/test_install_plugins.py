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
# Note: Due to the fact that it uses importlib, stevedore will see the first
# installed app immediately, but subsequent apps may not be seen until the 
# importlib cache is eventually updated or the interpreter reloaded.  
# importlib.invalidate_caches() does not appear to work as expected.


import importlib
import json
import os
from pathlib import Path
import site
import sys
import time

from contextlib import contextmanager
from inspect import getsourcelines
from subprocess import check_call


import pytest

from cogent3.app import APP_ENTRY_POINT, get_app_manager, get_app
from cogent3.app.composable import define_app
from cogent3.app import typing as c3types

def is_package_installed(package_name: str) -> bool:
    """Check if a package is installed"""
    site_packages_dir = site.getsitepackages()[1]
    installed_packages = os.listdir(site_packages_dir)
    return any(package_name in pkg for pkg in installed_packages)

def is_package_imported(package_name: str) -> bool:
    """Check if a package is imported"""
    try:
        importlib.import_module(package_name)
        return True
    except ImportError as e:
        return False
    
def is_app_installed(module_name: str, app_name: str) -> bool:
    """Check if app is installed"""
    app_manager = get_app_manager(force=True)
    return any(ext.entry_point.value == f"{module_name}:{app_name}" for ext in app_manager)    

@contextmanager
def persist_for(seconds:int, operation_name:str="Operation"):
    """A context manager that yields until the operation is complete or the time limit is reached."""
    start_time = time.time()
    try:
        def perform_operation(operation):
            while time.time() - start_time < seconds and not operation():
                time.sleep(1)
        yield perform_operation
    finally:
        elapsed_time = time.time() - start_time
        if elapsed_time > seconds:
            raise TimeoutError(f"{operation_name} timed out after {seconds} seconds.")

def install_app(cls, temp_dir: Path, mod=None):
    """
    Installs a temporary app created from a class, by writing the source code
    into a new Python file in the temporary directory, creating a setup.py file
    that references the new module, and running pip install . the temporary directory.
    """
    # Get the source code of the function
    lines, _ = getsourcelines(cls)

    # Find the first line of the class definition
    class_line = next(i for i, line in enumerate(lines) if line.strip().startswith("class"))
    
    # Calculate the number of leading spaces in the first line of the class definition
    leading_spaces = len(lines[class_line]) - len(lines[class_line].lstrip())


    # Remove leading spaces from all lines
    source_code = "\n".join(line[leading_spaces:] for line in lines)
    if not mod:
        mod = f"mod_{cls.__name__}"
    
    # Write the source code into a new Python file in the temporary directory
    package_path = temp_dir / f"{mod}.py"
    with open(package_path, "w") as package_file:
        package_file.write("from cogent3.app.composable import define_app\n")
        if hasattr(cls, "_requires_imports"):
            for import_statement in cls._requires_imports():
                package_file.write(f"{import_statement}\n")
        package_file.write(source_code)

    # Create a setup.py file in the temporary directory that references the new module
    setup_path = temp_dir / "setup.py"
    setup_contents = f"""
from setuptools import setup

setup(
name='{mod}',
version='0.1',
py_modules=['{mod}'],
entry_points={{
    "{APP_ENTRY_POINT}": [
        "{cls.__name__} = {mod}:{cls.__name__}",
    ],
}},    
)
"""
    with open(setup_path, "w") as setup_file:
        setup_file.write(setup_contents)

    # Use subprocess to run pip install . in the temporary directory
    check_call([sys.executable, "-m", "pip", "install", "."], cwd=temp_dir)

    with persist_for(seconds=60, operation_name=f"Installing package {mod}") as repeat_until:
        repeat_until(lambda: is_package_installed(package_name=mod))

    with persist_for(seconds=60, operation_name=f"Importing package {mod}") as repeat_until:
        repeat_until(lambda: is_package_imported(package_name=mod))

    with persist_for(seconds=60, operation_name=f"Importing package {mod}") as repeat_until:
        repeat_until(lambda: is_app_installed(module_name=mod,app_name=cls.__name__))

def uninstall_app(cls, mod=None):
    """
    Uninstalls a temporary app created from a class by running pip uninstall -y in the temporary directory.
    """
    if not mod:
        mod = f"mod_{cls.__name__}"
    package_name = mod.replace("_", "-")  # underscores in package names are replaced with hyphens
    check_call([sys.executable, "-m", "pip", "uninstall", "-y", package_name])
    _ = get_app_manager(force=True) # force the app cache to reload

@contextmanager
def load_app(app_class, tmp_path, module_name: str = None):
    """
    A context manager that creates a temporary app from a class, in a passed in path, 
    installs it, yields for tests then uninstalls the app.

    Usage
    -----
    with temp_app(cls,tmp_path):
        myapp = get_app(cls.__name__)
        # test myapp ... 
    """
    if module_name is None:
        module_name = f"mod_{app_class.__name__}"
    try:
        try:
            install_app(app_class, temp_dir=tmp_path, mod=module_name)
        except TimeoutError as e:
            pytest.fail(e.args[0])
        except Exception as e:
            pytest.fail(f"App installation failed: {e}")
        yield
    finally:
        uninstall_app(app_class, mod=module_name)

def test_install_temp_app(tmp_path: Path):
    @define_app
    class MyAppClass:
        def main(self, data: c3types.SerialisableType) -> str:
            return json.dumps(data)

        @classmethod
        def _requires_imports(cls):
            return ["from cogent3.app import typing as c3types", "import json"]

    with load_app(MyAppClass, tmp_path):
        myapp = get_app("MyAppClass")
        data = {"key": "value", "number": 42, "list": [1, 2, 3], "boolean": True}
        assert (
            myapp(data)
            == '{"key": "value", "number": 42, "list": [1, 2, 3], "boolean": true}'
        )