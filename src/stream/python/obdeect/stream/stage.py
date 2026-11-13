"""Write a pipeline stage in Python: ``transform`` rewrites the photons of a stream.

>>> def func(telescope_id, columns, metadata):
...     return {**columns, "w": columns["w"] * 0.9}
>>> transform("in.ps", "out.ps", func, stage_id=0x444C0001, frame_out="camera")
"""

from __future__ import annotations

import os
from collections.abc import Callable, Mapping

import numpy as np
import obdeect.stream as obdeect_stream
from obdeect.stream.stream import METADATA_RECORDS, Reader, Writer, _frame, reject_same_file

EFFICIENCY_STATUS_STATES = ("applied", "envelope", "pending")
_EFFICIENCY_STATUS_FIELDS = ("state", "source", "wavelength_nm", "value")


def efficiency_status_record(
    efficiency_status: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    """Build an ``efficiency_status`` record from a mapping of factor names to fields.

    ``state`` is ``applied``, ``envelope`` or ``pending``; ``source`` is the provenance string;
    ``wavelength_nm`` and ``value`` (equal-length float arrays) are optional.
    """
    metadata: dict[str, object] = {"factors": ",".join(efficiency_status)}
    for name, entry in efficiency_status.items():
        if not name or "," in name or "." in name:
            raise ValueError(f"invalid efficiency status factor name {name!r}")
        unknown = set(entry) - set(_EFFICIENCY_STATUS_FIELDS)
        if unknown:
            raise ValueError(f"efficiency status factor '{name}': unknown fields {sorted(unknown)}")
        if entry.get("state") not in EFFICIENCY_STATUS_STATES:
            raise ValueError(
                f"efficiency status factor '{name}': "
                f"state must be one of {EFFICIENCY_STATUS_STATES}"
            )
        if not isinstance(entry.get("source"), str):
            raise ValueError(f"efficiency status factor '{name}': source must be a string")
        if ("wavelength_nm" in entry) != ("value" in entry):
            raise ValueError(
                f"efficiency status factor '{name}': wavelength_nm and value go together"
            )
        if "value" in entry and np.size(entry["value"]) != np.size(entry["wavelength_nm"]):
            raise ValueError(
                f"efficiency status factor '{name}': wavelength_nm and value differ in length"
            )
        for field in _EFFICIENCY_STATUS_FIELDS:
            if field in entry:
                metadata[f"{name}.{field}"] = entry[field]
    return metadata


def _applied(metadata: Mapping[str, object]) -> set[str]:
    names = str(metadata.get("factors", "")).split(",")
    return {n for n in names if metadata.get(f"{n}.state") == "applied"}


def transform(
    input: str | os.PathLike,  # noqa: A002
    output: str | os.PathLike,
    func: Callable[[int, dict[str, np.ndarray], dict[str, object]], Mapping[str, np.ndarray]],
    *,
    stage_id: int,
    frame_out: int | str,
    efficiency_status: Mapping[str, Mapping[str, object]] | None = None,
) -> None:
    """Copy a stream, replacing the columns of every photon chunk with ``func``'s result.

    ``func(telescope_id, columns, metadata)`` gets the chunk's float32 columns (name to numpy
    array) and returns a mapping of the same kind for the output chunk: float32 arrays of
    equal length (at most ``obdeect.stream.max_chunk_rows``, may differ from the input;
    ``{}`` writes an empty chunk). ``metadata`` is a dict with the latest ``run``, ``event`` and
    ``telescope``
    record metadata (empty dicts before the first), the chunk's ``frame`` and ``first_index``,
    and ``stage_id``. Output indices are contiguous and the output uses ``frame_out``.
    ``telescope_end.photons`` is updated to the sum of output weights.

    All other records known to the reader are copied unchanged (records of unknown type are
    dropped, as the reader skips them), except that the ``frame`` of ``telescope`` records
    becomes ``frame_out``. ``efficiency_status`` (see ``efficiency_status_record``) is written
    as one record per run, after existing efficiency status records and before the first event.
    A factor ``applied`` in any input record, including one after a ``telescope`` record,
    cannot be applied again (ValueError). The writer's
    ``producer`` is ``"obdeect.stream.stage 0x<stage_id>"``.
    """
    if isinstance(stage_id, bool) or not isinstance(stage_id, int) or not 0 <= stage_id < 2**32:
        raise ValueError("stage_id must be an integer in [0, 2**32)")
    frame_code = _frame(frame_out)
    if frame_code not in obdeect_stream.FRAMES.values():
        raise ValueError(f"unknown frame {frame_out!r}")
    entry = efficiency_status_record(efficiency_status) if efficiency_status else None
    context: dict[str, dict[str, object]] = {"run": {}, "event": {}, "telescope": {}}
    entry_applied = _applied(entry) if entry else set()
    efficiency_status_pending = False
    output_index: dict[int, int] = {}  # per telescope; continues across array reuses
    output_weight = 0.0

    reject_same_file(input, output)
    with Writer(os.fspath(output), f"obdeect.stream.stage {stage_id:#x}") as writer:
        for record in Reader(os.fspath(input)):
            kind = record["type"]
            if kind == "stream":
                continue
            if efficiency_status_pending and kind in ("event", "run_end"):
                writer.write("efficiency_status", entry)
                efficiency_status_pending = False
            if kind in METADATA_RECORDS:
                metadata = dict(record["metadata"])
                if kind in context:
                    context[kind] = metadata
                if kind == "event":
                    output_index = {}
                elif kind == "telescope":
                    output_weight = 0.0
                elif kind == "telescope_end":
                    metadata["photons"] = output_weight
                if kind == "run":
                    efficiency_status_pending = entry is not None
                elif kind == "efficiency_status":
                    repeated = entry_applied & _applied(metadata)
                    if repeated:
                        raise ValueError(f"factor(s) already applied upstream: {sorted(repeated)}")
                elif kind == "telescope" and "frame" in metadata:
                    metadata["frame"] = frame_code
                writer.write(kind, metadata)
            elif kind == "blob":
                writer.write_blob(record["eventio_type"], record["eventio_id"], record["data"])
            elif kind == "photons":
                info = {
                    **context,
                    "frame": record["frame"],
                    "first_index": record["first_index"],
                    "stage_id": stage_id,
                }
                result = func(record["telescope_id"], record["columns"], info)
                telescope = record["telescope_id"]
                writer.write_photons(telescope, frame_code, output_index.get(telescope, 0), result)
                output_index[telescope] = output_index.get(telescope, 0) + len(result.get("x", ()))
                output_weight += float(np.sum(result.get("w", ()), dtype=np.float64))
