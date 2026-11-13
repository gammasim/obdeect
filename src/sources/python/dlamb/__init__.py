"""dlamb: artificial light sources for Cherenkov-telescope simulations."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("dlamb")
except PackageNotFoundError:
    __version__ = "0.0.0"
