"""gstrickd: photon-stream format and I/O for Cherenkov-telescope simulation chains."""

from importlib.metadata import PackageNotFoundError, version

from gstrickd.stream import COLUMNS, FRAMES, METADATA_RECORDS, REQUIRED_COLUMNS, Reader, Writer

try:
    __version__ = version("gstrickd")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = [
    "COLUMNS",
    "FRAMES",
    "METADATA_RECORDS",
    "REQUIRED_COLUMNS",
    "Reader",
    "Writer",
    "__version__",
]
