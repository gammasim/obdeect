"""dlufd: atmospheric propagation of Cherenkov and artificial-source photons."""

from importlib.metadata import PackageNotFoundError, version

from dlufd.profile import (
    AtmosphereProfile,
    RefractionTables,
    air_delay_ns,
    propagation_time_ns,
    refract,
)

try:
    __version__ = version("dlufd")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = [
    "AtmosphereProfile",
    "RefractionTables",
    "__version__",
    "air_delay_ns",
    "propagation_time_ns",
    "refract",
]
