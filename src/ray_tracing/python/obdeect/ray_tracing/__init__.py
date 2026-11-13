"""obdeect: optical ray tracing for imaging atmospheric Cherenkov telescopes."""

from obdeect.stream._package import lazy_exports, package_version

__version__ = package_version("obdeect-ray-tracing")

# Imported on first use, so that reading the version needs no NumPy.
_EXPORTS = {
    "DETECTED": "obdeect.ray_tracing.model",
    "STATUS_NAMES": "obdeect.ray_tracing.model",
    "OpticalModel": "obdeect.ray_tracing.model",
    "ThreadPool": "obdeect.ray_tracing.model",
    "trace_many": "obdeect.ray_tracing.model",
    "ModelError": "obdeect.ray_tracing.formats",
}
__all__ = list(_EXPORTS)

__getattr__ = lazy_exports("obdeect", _EXPORTS)
