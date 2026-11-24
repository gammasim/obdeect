"""The package and its native extension import; the public names resolve; imports are lazy."""

import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version

import obdeect.propagation as obdeect_propagation
import pytest
from obdeect.propagation import _core


def test_version():
    assert isinstance(obdeect_propagation.__version__, str)
    assert len(_core.core_version()) == 3


def test_version_is_the_installed_distribution_version():
    try:
        installed = version("obdeect-propagation")
    except PackageNotFoundError:
        pytest.fail("obdeect-propagation is not installed")
    assert obdeect_propagation.__version__ == installed


def test_every_public_name_resolves():
    assert "propagate" in obdeect_propagation.__all__
    for name in obdeect_propagation.__all__:
        assert getattr(obdeect_propagation, name) is not None, name


def test_unknown_attribute_is_an_error():
    with pytest.raises(AttributeError, match="no attribute 'not_an_export'"):
        getattr(obdeect_propagation, "not_an_export")


def test_version_import_does_not_load_the_extension():
    # A fresh interpreter: reading the version must not import NumPy or the native core.
    code = (
        "import sys, obdeect.propagation as p; "
        "p.__version__; "
        "assert 'obdeect.propagation._core' not in sys.modules, 'core loaded eagerly'; "
        "assert 'numpy' not in sys.modules, 'numpy loaded eagerly'"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
