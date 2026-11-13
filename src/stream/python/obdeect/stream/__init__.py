"""gstrickd: photon-stream format and I/O for Cherenkov-telescope simulation chains."""

import importlib
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("obdeect-stream")
except PackageNotFoundError:
    __version__ = "0.0.0"

# Imported on first use, so that reading the version needs neither NumPy nor the extension.
# The format constants are read from the extension (gstrickd._core).
_EXPORTS = {
    "COLUMNS": "obdeect.stream.stream",
    "column_names": "obdeect.stream._core",
    "corsika_to_stream": "obdeect.stream.convert",
    "FRAMES": "obdeect.stream.stream",
    "frame_ids": "obdeect.stream._core",
    "format_major": "obdeect.stream._core",
    "format_minor": "obdeect.stream._core",
    "max_chunk_rows": "obdeect.stream._core",
    "METADATA_RECORDS": "obdeect.stream.stream",
    "REQUIRED_COLUMNS": "obdeect.stream.stream",
    "stream_to_grisu": "obdeect.stream.convert",
    "stream_to_simtel": "obdeect.stream.convert",
    "Reader": "obdeect.stream.stream",
    "Writer": "obdeect.stream.stream",
}
__all__ = [*_EXPORTS, "__version__"]


def __getattr__(name: str):
    if name in _EXPORTS:
        return getattr(importlib.import_module(_EXPORTS[name]), name)
    raise AttributeError(f"module 'gstrickd' has no attribute {name!r}")
