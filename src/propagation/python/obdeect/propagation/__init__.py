"""obdeect propagation: atmospheric propagation of Cherenkov and artificial-source photons."""

from obdeect.stream._package import lazy_exports, package_version

__version__ = package_version("obdeect-propagation")

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

__getattr__ = lazy_exports("obdeect.propagation", _EXPORTS)
