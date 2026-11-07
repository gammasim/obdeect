"""dlufd: atmospheric propagation of Cherenkov and artificial-source photons."""

import importlib
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("dlufd")
except PackageNotFoundError:
    __version__ = "0.0.0"

# Imported on first use, so that reading the version needs neither NumPy nor the extension.
_EXPORTS = {
    "TransmissionTable": "dlufd.extinction",
    "apply_extinction": "dlufd.extinction",
    "AtmosphereProfile": "dlufd.profile",
    "RefractionTables": "dlufd.profile",
    "air_delay_ns": "dlufd.profile",
    "propagation_time_ns": "dlufd.profile",
    "refract": "dlufd.profile",
    "propagate": "dlufd.propagation",
}
__all__ = [*_EXPORTS, "__version__"]


def __getattr__(name: str):
    if name in _EXPORTS:
        return getattr(importlib.import_module(_EXPORTS[name]), name)
    raise AttributeError(f"module 'dlufd' has no attribute {name!r}")
