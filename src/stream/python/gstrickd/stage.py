"""Write a pipeline stage in Python: ``transform`` rewrites the photons of a stream.

>>> def func(telescope_id, columns, metadata):
...     return {**columns, "w": columns["w"] * 0.9}
>>> transform("in.ps", "out.ps", func, stage_id=0x444C0001, frame_out="camera")
"""

from __future__ import annotations

import os
from collections.abc import Callable, Mapping

import numpy as np

import gstrickd
from gstrickd.stream import METADATA_RECORDS, Reader, Writer, _frame

LEDGER_STATES = ("applied", "envelope", "pending")
_LEDGER_FIELDS = ("state", "source", "wavelength_nm", "value")


def ledger_record(ledger: Mapping[str, Mapping[str, object]]) -> dict[str, object]:
    """Metadata of a ``ledger`` record from ``{factor: {state, source, wavelength_nm, value}}``.

    ``state`` is ``applied``, ``envelope`` or ``pending``; ``source`` is the provenance string;
    ``wavelength_nm`` and ``value`` (equal-length float arrays) are optional.
    """
    metadata: dict[str, object] = {"factors": ",".join(ledger)}
    for name, entry in ledger.items():
        if not name or "," in name or "." in name:
            raise ValueError(f"invalid ledger factor name {name!r}")
        unknown = set(entry) - set(_LEDGER_FIELDS)
        if unknown:
            raise ValueError(f"ledger factor '{name}': unknown fields {sorted(unknown)}")
        if entry.get("state") not in LEDGER_STATES:
            raise ValueError(f"ledger factor '{name}': state must be one of {LEDGER_STATES}")
        if not isinstance(entry.get("source"), str):
            raise ValueError(f"ledger factor '{name}': source must be a string")
        if ("wavelength_nm" in entry) != ("value" in entry):
            raise ValueError(f"ledger factor '{name}': wavelength_nm and value go together")
        if "value" in entry and np.size(entry["value"]) != np.size(entry["wavelength_nm"]):
            raise ValueError(f"ledger factor '{name}': wavelength_nm and value differ in length")
        for field in _LEDGER_FIELDS:
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
    ledger: Mapping[str, Mapping[str, object]] | None = None,
) -> None:
    """Copy a stream, replacing the columns of every photon chunk with ``func``'s result.

    ``func(telescope_id, columns, metadata)`` gets the chunk's float32 columns (name to numpy
    array) and returns a mapping of the same kind for the output chunk: float32 arrays of
    equal length (at most ``gstrickd.max_chunk_rows``, may differ from the input; ``{}`` writes
    an empty chunk). ``metadata`` is a dict with the latest ``run``, ``event`` and ``telescope``
    record metadata (empty dicts before the first), the chunk's ``frame`` and ``first_index``,
    and ``stage_id``. The output chunk keeps ``first_index`` and gets ``frame_out``.

    All other records known to the reader are copied unchanged (records of unknown type are
    dropped, as the reader skips them), except that the ``frame`` of ``telescope`` records
    becomes ``frame_out``. ``ledger`` (see ``ledger_record``) is written as one ``ledger``
    record per run, after the ledger records already there and before the first event. A
    factor ``applied`` in any input ledger, also one that follows a ``telescope`` record,
    cannot be applied again (ValueError). The writer's
    ``producer`` is ``"gstrickd.stage 0x<stage_id>"``.
    """
    if isinstance(stage_id, bool) or not isinstance(stage_id, int) or not 0 <= stage_id < 2**32:
        raise ValueError("stage_id must be an integer in [0, 2**32)")
    frame_code = _frame(frame_out)
    if frame_code not in gstrickd.FRAMES.values():
        raise ValueError(f"unknown frame {frame_out!r}")
    entry = ledger_record(ledger) if ledger else None
    context: dict[str, dict[str, object]] = {"run": {}, "event": {}, "telescope": {}}
    entry_applied = _applied(entry) if entry else set()
    ledger_pending = False

    with Writer(os.fspath(output), f"gstrickd.stage {stage_id:#x}") as writer:
        for record in Reader(os.fspath(input)):
            kind = record["type"]
            if kind == "stream":
                continue
            if ledger_pending and kind in ("event", "run_end"):
                writer.write("ledger", entry)
                ledger_pending = False
            if kind in METADATA_RECORDS:
                metadata = dict(record["metadata"])
                if kind in context:
                    context[kind] = metadata
                if kind == "run":
                    ledger_pending = entry is not None
                elif kind == "ledger":
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
                writer.write_photons(
                    record["telescope_id"], frame_code, record["first_index"], result
                )
