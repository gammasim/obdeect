"""gstrickd: photon-stream format and I/O for Cherenkov-telescope simulation chains."""

import importlib
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("gstrickd")
except PackageNotFoundError:
    __version__ = "0.0.0"

# Imported on first use, so that reading the version needs neither NumPy nor the extension.
_EXPORTS = {
    "COLUMNS": "gstrickd.stream",
    "FRAMES": "gstrickd.stream",
    "METADATA_RECORDS": "gstrickd.stream",
    "REQUIRED_COLUMNS": "gstrickd.stream",
    "Reader": "gstrickd.stream",
    "Writer": "gstrickd.stream",
}
__all__ = [*_EXPORTS, "__version__"]


def __getattr__(name: str):
    if name in _EXPORTS:
        return getattr(importlib.import_module(_EXPORTS[name]), name)
    raise AttributeError(f"module 'gstrickd' has no attribute {name!r}")
