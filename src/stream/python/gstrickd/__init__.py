"""gstrickd: photon-stream format and I/O for Cherenkov-telescope simulation chains."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("gstrickd")
except PackageNotFoundError:
    __version__ = "0.0.0"
