"""obdeect propagation: atmospheric propagation of Cherenkov and artificial-source photons."""

import importlib
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("obdeect-propagation")
except PackageNotFoundError:
    __version__ = "0.0.0"

# Imported on first use, so that reading the version needs neither NumPy nor the extension.
_EXPORTS = {
    "TransmissionTable": "obdeect.propagation.extinction",
    "apply_extinction": "obdeect.propagation.extinction",
    "AtmosphereProfile": "obdeect.propagation.profile",
    "ScatteringTable": "obdeect.propagation.scattering",
    "scatter": "obdeect.propagation.scattering",
    "RefractionTables": "obdeect.propagation.profile",
    "air_delay_ns": "obdeect.propagation.profile",
    "propagation_time_ns": "obdeect.propagation.profile",
    "refract": "obdeect.propagation.profile",
    "propagate": "obdeect.propagation.propagation",
    "convert_table": "obdeect.propagation.convert",
}
__all__ = [*_EXPORTS, "__version__"]


def __getattr__(name: str):
    if name in _EXPORTS:
        return getattr(importlib.import_module(_EXPORTS[name]), name)
    raise AttributeError(f"module 'obdeect.propagation' has no attribute {name!r}")
