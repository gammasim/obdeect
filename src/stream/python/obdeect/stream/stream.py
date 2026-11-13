"""Photon-stream API: Reader and Writer with Python-side conveniences."""

from __future__ import annotations

import os
from collections.abc import Mapping

import numpy as np
from obdeect.stream import _core

Reader = _core.Reader

FRAMES = dict(_core.frame_ids)
COLUMNS = tuple(_core.column_names())
REQUIRED_COLUMNS = ("x", "y", "cx", "cy", "t", "w")
METADATA_RECORDS = (
    "run",
    "event",
    "telescope",
    "telescope_end",
    "event_end",
    "run_end",
    "efficiency_status",
)


def _frame(frame: int | str) -> int:
    if isinstance(frame, str):
        if frame not in FRAMES:
            raise ValueError(f"unknown frame '{frame}', expected one of {sorted(FRAMES)}")
        return FRAMES[frame]
    if isinstance(frame, (bool, np.bool_)) or not isinstance(frame, (int, np.integer)):
        raise TypeError(f"frame must be a name or an integer, not {frame!r}")
    return int(frame)


def _metadata_value(value):
    if isinstance(value, (bool, np.bool_, int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value)
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple, np.ndarray)):
        return np.asarray(value, dtype=np.float64).ravel().tolist()
    raise TypeError(f"unsupported metadata value of type {type(value).__name__}")


def reject_same_file(input: object, output: object) -> None:
    """Raise ValueError if output is the same file as input (also via links)."""
    source, target = os.fspath(input), os.fspath(output)
    if source == "-" or target == "-":
        return
    try:
        same = os.path.samefile(source, target)
    except OSError:
        return  # one of them does not exist yet
    if same:
        raise ValueError(f"input and output are the same file: {target}")


class Writer(_core.Writer):
    """Writes a photon stream. As a context manager it writes stream_end on a normal exit."""

    def __enter__(self) -> Writer:
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        if exc_type is None:
            self.close()
        else:
            self.abort()  # no stream_end: the partial output is rejected as truncated

    def write(self, record_type: str, metadata: Mapping[str, object]) -> None:
        """Write a metadata record, e.g. write("run", {"run_number": 1})."""
        super().write(record_type, {k: _metadata_value(v) for k, v in metadata.items()})

    def write_photons(
        self,
        telescope_id: int,
        frame: int | str,
        first_index: int,
        columns: Mapping[str, np.ndarray],
    ) -> None:
        """Write one photon chunk; columns map names (see COLUMNS) to 1-D float32 arrays."""
        arrays = {}
        for name, values in columns.items():
            array = np.asarray(values)
            if array.dtype != np.float32:
                raise TypeError(f"column '{name}' must be float32, got {array.dtype}")
            arrays[name] = np.ascontiguousarray(array)
        super().write_photons(telescope_id, _frame(frame), first_index, arrays)


__all__ = ["COLUMNS", "FRAMES", "METADATA_RECORDS", "REQUIRED_COLUMNS", "Reader", "Writer"]
