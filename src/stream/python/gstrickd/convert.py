"""Conversions at the boundary: CORSIKA IACT output to a stream, and a stream to sim_telarray."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping

from gstrickd import _core

MIRROR_CLASSES = {1: "single mirror (Bypass_Optics=1)", 2: "dual mirror (Bypass_Optics=2)"}


def corsika_to_stream(input_path: str | os.PathLike, output_path: str | os.PathLike) -> None:
    """Convert CORSIKA IACT EventIO (file or ``"-"`` for stdin) to a photon stream file."""
    _core.corsika_to_stream(os.fspath(input_path), os.fspath(output_path))


def _camera_entries(camera_map) -> list[tuple[int, int, float]]:
    if camera_map is None:
        return []
    if isinstance(camera_map, (str, os.PathLike)):
        with open(camera_map, encoding="utf-8") as file:
            return list(_core.parse_camera_map(file.read()))
    if not isinstance(camera_map, Mapping):
        raise TypeError(
            "camera_map must be None, a path, or a mapping of telescope id to "
            "(mirror_class, rotation_deg)"
        )
    entries = []
    for telescope_id, value in camera_map.items():
        mirror_class, rotation_deg = value
        entries.append((int(telescope_id), int(mirror_class), float(rotation_deg)))
    return entries


def stream_to_simtel(
    input_path: str | os.PathLike,
    output_path: str | os.PathLike,
    camera_map: Mapping[int, tuple[int, float]] | str | os.PathLike | None = None,
    allow_pending: Iterable[str] = (),
) -> None:
    """Convert a photon stream file to sim_telarray input (EventIO).

    ``camera_map`` lists the camera-frame telescopes: a mapping of telescope id (1-based) to
    ``(mirror_class, camera_rotation_deg)``, or the path of a text file with one line
    ``telescope_id mirror_class camera_rotation_deg`` per telescope. Telescopes not listed are
    ground-frame telescopes (Bypass_Optics=0 in sim_telarray).

    The efficiency ledger must list only factors that sim_telarray applies itself in the
    telescope's frame as ``pending`` or ``envelope`` (docs/usage.md). Other factors raise
    ``RuntimeError``; ``allow_pending`` names factors to accept as they are.
    """
    if isinstance(allow_pending, str):
        raise TypeError("allow_pending must be a sequence of factor names, not a string")
    _core.stream_to_simtel(
        os.fspath(input_path),
        os.fspath(output_path),
        _camera_entries(camera_map),
        [str(name) for name in allow_pending],
    )
