"""obdeect stream: photon-stream format and I/O for Cherenkov-telescope simulation chains."""

from obdeect.stream._package import lazy_exports, package_version

__version__ = package_version("obdeect-stream")

# Imported on first use, so that reading the version needs neither NumPy nor the extension.
# The format constants are read from the extension (obdeect.stream._core).
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


__getattr__ = lazy_exports("obdeect.stream", _EXPORTS)
