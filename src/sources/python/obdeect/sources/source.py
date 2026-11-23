"""Photons of single emissions: the emission of a run, without telescopes or events.

An emission is a dictionary with a ``kind`` (point, line, plane_wave or sky_background) and
the keys of :mod:`obdeect.sources.config`. The functions here accept a sparse dictionary and
normalize it first; ``normalize_emission`` returns the normalized form.

.. code-block:: python

    flasher = {
        "kind": "point",
        "photons": 1e6,  # emitted into 4 pi
        "position_m": [0.0, 0.0, 1.5],
        "axis": [0.0, 0.0, -1.0],
        "angular": {"kind": "cone", "angle_deg": 11.0},
        "spectrum": {"kind": "fixed", "wavelength_nm": 385.0},
        "pulse": {"kind": "tophat", "width_ns": 2.5},
        "target": {"kind": "disc", "centre_m": [0.0, 0.0, 0.0],
                   "normal": [0.0, 0.0, 1.0], "radius_m": 1.2},
    }
    photons = emit(flasher, seed=3)  # {"x": float32 array, ..., "w": ..., "wavelength": ...}

An ``auto`` target needs the telescopes of a run and is rejected here.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import Any

import numpy as np
from obdeect.sources import _core, config

# Copies of obdeect stream constants (no runtime dependency on obdeect stream). The tests in
# python/tests/test_stream_contract.py compare them with the compiled obdeect stream headers.

#: Largest number of photons in one chunk (photon stream chunk limit).
CHUNK_ROWS = 1 << 16

#: Frames of the photon stream (G is not a source frame); a run sets its frame.
FRAMES = {"sky": 0, "telescope": 2, "camera": 3}

#: Columns of the emitted photons, in obdeect stream column names.
COLUMNS = ("x", "y", "z", "cx", "cy", "cz", "t", "w", "wavelength")

_EVENT_PATH = "events[0].emissions[0]"


def normalize_emission(spec: Mapping[str, Any]) -> dict[str, Any]:
    """Normalized emission (defaults filled in); errors name the key within the emission."""
    if not isinstance(spec, Mapping):
        raise ValueError(f"emission: expected a dictionary, got {type(spec).__name__}")
    wrapper = {
        "run": {"observation_level_m": 0.0},
        "telescopes": [{"id": 1, "radius_m": 1.0}],
        "events": [{"emissions": [dict(spec)]}],
    }
    try:
        return config.normalize(wrapper)["events"][0]["emissions"][0]
    except ValueError as error:
        message = str(error).replace(_EVENT_PATH, "emission")
        raise ValueError(message) from None


def _check_rows(rows: int) -> None:
    if not 0 < rows <= CHUNK_ROWS:
        raise ValueError(f"rows must be in 1..{CHUNK_ROWS}")


def expected_photons(spec: Mapping[str, Any]) -> float:
    """Mean number of candidate photons of an emission (before line photons miss a target)."""
    return _core.expected_count(normalize_emission(spec))


def expected_weight(spec: Mapping[str, Any]) -> float:
    """Expected sum of the photon weights of an emission (its bunch size included)."""
    return _core.expected_weight(normalize_emission(spec))


def photon_count(spec: Mapping[str, Any], *, seed: int, event: int = 0, telescope: int = 0) -> int:
    """Number of candidate photons of an emission in one (event, telescope), keyed on the seed."""
    return _core.photon_count(normalize_emission(spec), seed, event, telescope)


def emit_chunk(
    spec: Mapping[str, Any],
    *,
    seed: int,
    event: int = 0,
    telescope: int = 0,
    first_index: int = 0,
    rows: int = CHUNK_ROWS,
) -> dict[str, np.ndarray]:
    """Candidates first_index .. first_index + rows - 1 (at most 2^16) as float32 columns.

    Line photons that miss their target, and photons of a ``to_plane_z0`` emission that move away
    from the plane z = 0, are dropped, so the chunk may hold fewer rows.
    """
    _check_rows(rows)
    if first_index < 0:
        raise ValueError("first_index must be >= 0")
    return _core.emit_chunk(normalize_emission(spec), seed, event, telescope, first_index, rows)


def iter_chunks(
    spec: Mapping[str, Any], *, seed: int, event: int = 0, telescope: int = 0
) -> Iterator[dict[str, np.ndarray]]:
    """All photons of an emission in chunks of at most 2^16 candidates (bounded memory)."""
    emission = normalize_emission(spec)
    count = _core.photon_count(emission, seed, event, telescope)
    for first in range(0, count, CHUNK_ROWS):
        rows = min(CHUNK_ROWS, count - first)
        yield _core.emit_chunk(emission, seed, event, telescope, first, rows)


def emit(
    spec: Mapping[str, Any], *, seed: int, event: int = 0, telescope: int = 0
) -> dict[str, np.ndarray]:
    """All photons of an emission as float32 columns (obdeect stream names; see ``COLUMNS``).

    The arrays hold every photon in memory; use ``iter_chunks`` for large counts.
    """
    chunks = list(iter_chunks(spec, seed=seed, event=event, telescope=telescope))
    if not chunks:
        return {name: np.empty(0, dtype=np.float32) for name in COLUMNS}
    return {name: np.concatenate([chunk[name] for chunk in chunks]) for name in COLUMNS}
