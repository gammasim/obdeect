"""Helpers shared by the obdeect component packages: installed version and lazy exports."""

import importlib
from collections.abc import Callable, Mapping
from importlib.metadata import PackageNotFoundError, version


def package_version(distribution: str) -> str:
    """Installed version of a distribution; ``0.0.0`` if it is not installed."""
    try:
        return version(distribution)
    except PackageNotFoundError:
        return "0.0.0"


def lazy_exports(module_name: str, exports: Mapping[str, str]) -> Callable[[str], object]:
    """Module ``__getattr__`` that imports each export from its module on first use.

    ``exports`` maps an attribute name to the dotted name of the module that defines it.
    """

    def __getattr__(name: str):
        if name in exports:
            return getattr(importlib.import_module(exports[name]), name)
        raise AttributeError(f"module {module_name!r} has no attribute {name!r}")

    return __getattr__
