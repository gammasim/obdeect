"""gstrickd: photon-stream format and I/O for Cherenkov-telescope simulation chains."""

import importlib
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("gstrickd")
except PackageNotFoundError:
    __version__ = "0.0.0"

# Imported on first use, so that reading the version needs neither NumPy nor the extension.
# The format constants are read from the extension (gstrickd._core).
_EXPORTS = {
    "COLUMNS": "gstrickd.stream",
    "column_names": "gstrickd._core",
    "corsika_to_stream": "gstrickd.convert",
    "FRAMES": "gstrickd.stream",
    "frame_ids": "gstrickd._core",
    "format_major": "gstrickd._core",
    "format_minor": "gstrickd._core",
    "max_chunk_rows": "gstrickd._core",
    "METADATA_RECORDS": "gstrickd.stream",
    "REQUIRED_COLUMNS": "gstrickd.stream",
    "stream_to_grisu": "gstrickd.convert",
    "stream_to_simtel": "gstrickd.convert",
    "Reader": "gstrickd.stream",
    "Writer": "gstrickd.stream",
}
__all__ = [*_EXPORTS, "__version__"]


def __getattr__(name: str):
    if name in _EXPORTS:
        return getattr(importlib.import_module(_EXPORTS[name]), name)
    raise AttributeError(f"module 'gstrickd' has no attribute {name!r}")
