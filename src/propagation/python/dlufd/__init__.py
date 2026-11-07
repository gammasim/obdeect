"""dlufd: atmospheric propagation of Cherenkov and artificial-source photons."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("dlufd")
except PackageNotFoundError:
    __version__ = "0.0.0"
