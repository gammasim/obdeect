"""Write a pipeline stage in Python: ``transform`` rewrites the photons of a stream.

>>> def func(telescope_id, columns, metadata):
...     return {**columns, "w": columns["w"] * 0.9}
>>> transform("in.pst", "out.pst", func, stage_id=0x444C0001, frame_out="camera")
"""

from __future__ import annotations

import os
from collections.abc import Callable, Mapping

import numpy as np
import obdeect.stream as obdeect_stream
from obdeect.stream.stream import (
    METADATA_RECORDS,
    Reader,
    Writer,
    _frame,
    reject_same_file,
)

EFFICIENCY_STATUS_STATES = ("applied", "envelope", "pending")
_EFFICIENCY_STATUS_FIELDS = ("state", "source", "wavelength_nm", "value")
_STATE_RANK = {"pending": 0, "envelope": 1, "applied": 2}


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


def _states(metadata: Mapping[str, object]) -> dict[str, str]:
    """Factor -> state of an efficiency_status record (the reader has validated the record)."""
    names = [n.strip() for n in str(metadata.get("factors", "")).split(",") if n.strip()]
    return {n: str(metadata[f"{n}.state"]) for n in names}


def _applied(metadata: Mapping[str, object]) -> set[str]:
    return {name for name, state in _states(metadata).items() if state == "applied"}


def _check_advances(run_states: Mapping[str, str], entry: Mapping[str, object]) -> None:
    """Raise ValueError unless each factor of ``entry`` advances from its input run state.

    A factor advances pending -> envelope -> applied (docs/stream/efficiency-status.md), so a
    factor already at the same or a later state in the input cannot be written again.
    """
    for name, state in _states(entry).items():
        previous = run_states.get(name)
        if previous is not None and _STATE_RANK[state] <= _STATE_RANK[previous]:
            raise ValueError(
                f"efficiency status factor '{name}' cannot go from {previous} to {state} "
                "(the input already has it in run scope)"
            )


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
    cannot be applied again, and a factor cannot advance to a state the input already has in
    run scope (ValueError). The writer's ``producer`` is ``"obdeect.stream.stage 0x<stage_id>"``.
    """
    if isinstance(stage_id, bool) or not isinstance(stage_id, int) or not 0 <= stage_id < 2**32:
        raise ValueError("stage_id must be an integer in [0, 2**32)")
    frame_code = _frame(frame_out)
    if frame_code not in obdeect_stream.FRAMES.values():
        raise ValueError(f"unknown frame {frame_out!r}")
    entry = efficiency_status_record(efficiency_status) if efficiency_status else None
    context: dict[str, dict[str, object]] = {"run": {}, "event": {}, "telescope": {}}
    entry_applied = _applied(entry) if entry else set()
    run_states: dict[str, str] = {}  # run-scope states of the input, before the first event
    seen_event = False
    efficiency_status_pending = False
    output_index: dict[int, int] = {}  # per telescope; continues across array reuses
    output_weight = 0.0

    reject_same_file(input, output)
    # The input is opened before the output is created (the writer truncates it), so a missing
    # input leaves an existing output alone.
    with (
        Reader(os.fspath(input)) as reader,
        Writer(os.fspath(output), f"obdeect.stream.stage {stage_id:#x}") as writer,
    ):
        for record in reader:
            kind = record["type"]
            if kind == "stream":
                continue
            if efficiency_status_pending and kind in ("event", "run_end"):
                _check_advances(run_states, entry)
                writer.write("efficiency_status", entry)
                efficiency_status_pending = False
            if kind in METADATA_RECORDS:
                metadata = dict(record["metadata"])
                if kind in context:
                    context[kind] = dict(metadata)  # a copy: the output frame is set below
                if kind == "event":
                    seen_event = True
                    output_index = {}
                elif kind == "telescope":
                    output_weight = 0.0
                elif kind == "telescope_end":
                    metadata["photons"] = output_weight
                if kind == "run":
                    run_states = {}
                    seen_event = False
                    efficiency_status_pending = entry is not None
                elif kind == "efficiency_status":
                    repeated = entry_applied & _applied(metadata)
                    if repeated:
                        raise ValueError(f"factor(s) already applied upstream: {sorted(repeated)}")
                    if not seen_event:
                        run_states.update(_states(metadata))
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
