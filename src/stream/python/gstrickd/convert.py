"""Conversions at the boundary: CORSIKA IACT output to a stream, and a stream to sim_telarray."""

from __future__ import annotations

import os
import warnings
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
    scaled_factors: Iterable[str] = (),
) -> None:
    """Convert a photon stream file to sim_telarray input (EventIO).

    ``camera_map`` lists the camera-frame telescopes: a mapping of telescope id (1-based) to
    ``(mirror_class, camera_rotation_deg)``, or the path of a text file with one line
    ``telescope_id mirror_class camera_rotation_deg`` per telescope. Telescopes not listed are
    ground-frame telescopes (Bypass_Optics=0 in sim_telarray).

    The efficiency ledger is checked per telescope against the frame of its photons
    (docs/format.md): a factor is valid as ``applied`` when sim_telarray does not apply it in
    that frame (e.g. ``atmosphere`` for camera-frame photons), and as ``pending`` when
    sim_telarray applies it (e.g. ``qe``). ``envelope`` is never valid for factors sim_telarray
    applies. With ground-frame photons ``atmosphere`` is valid in any state. Anything else, and
    unknown factor names, raise ``RuntimeError``; ``allow_pending`` names factors to accept as
    they are. ``scaled_factors`` names factors whose sim_telarray table holds f/f_env (e.g. a
    QE table divided by the QE envelope): they must be in state ``envelope``.
    """
    for name, value in (("allow_pending", allow_pending), ("scaled_factors", scaled_factors)):
        if isinstance(value, str):
            raise TypeError(f"{name} must be a sequence of factor names, not a string")
    _core.stream_to_simtel(
        os.fspath(input_path),
        os.fspath(output_path),
        _camera_entries(camera_map),
        [str(name) for name in allow_pending],
        [str(name) for name in scaled_factors],
    )


def stream_to_grisu(
    input_path: str | os.PathLike,
    output_path: str | os.PathLike,
    telescope_config: str | os.PathLike | None = None,
    qe: float = 1.0,
    expand_bunches: bool = False,
    event_details: bool = False,
) -> None:
    """Convert a ground-frame photon stream file to GrISU/GrOptics photon text input.

    ``telescope_config`` is a file with ``* TLLOC id x y z offset px py`` lines (x east, y
    north, metres): stream telescopes map to the GrISU id of the line within 0.5 m of their
    position; unmatched telescopes are skipped with a warning. Without it the GrISU id is the
    stream telescope id. ``qe`` is written as the global QE factor (``R`` line) and not
    applied. Photon bunches with weight above 1 raise ``RuntimeError`` unless
    ``expand_bunches`` writes ``round(w)`` copies. The wavelength must be assigned and the
    ledger factor ``atmosphere`` must be ``applied``; otherwise ``RuntimeError`` is raised.
    ``event_details`` adds a ``C`` line after each ``S`` line.
    """
    text = _core.stream_to_grisu(
        os.fspath(input_path),
        os.fspath(output_path),
        "" if telescope_config is None else os.fspath(telescope_config),
        float(qe),
        bool(expand_bunches),
        bool(event_details),
    )
    if text:
        warnings.warn(text.strip(), stacklevel=2)
