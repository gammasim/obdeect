"""Light sources from plain dictionaries, and emission of their photons.

A source is a dictionary with the keys below. Values are SI-like units (metres, ns,
nm, degrees); the frame and direction conventions are those of the photon stream
(gstrickd): ``cz < 0`` is downward, and the frames are ``sky``, ``telescope`` (T) and
``camera`` (C). Exactly one of ``emitter`` or ``star`` is given.

.. code-block:: python

    flasher = {
        "frame": "camera",
        "photons": 1e6,  # emitted into 4 pi
        "emitter": {
            "position_m": [0.0, 0.0, 1.5],
            "pointing": [0.0, 0.0, 1.0],
            "angular": {"kind": "filled_cone", "angle_deg": 11.0},
        },
        "spectrum": {"kind": "fixed", "wavelength_nm": 385.0},
        "pulse": {"kind": "tophat", "width_ns": 2.5, "delay_ns": 0.0},
        "target": {"kind": "disc", "centre_m": [0.0, 0.0, 0.0],
                   "normal": [0.0, 0.0, 1.0], "radius_m": 1.2},
    }
    photons = emit(flasher, seed=3)  # {"x": float32 array, ..., "w": ...}

Keys: ``frame``, ``photons``, ``emitter`` or ``star``, ``spectrum`` and ``pulse``
(default: fixed 400 nm and fixed time 0), ``target`` (optional, emitters only).
Unknown keys and kinds are errors.
"""

from __future__ import annotations

import math
from collections.abc import Iterator, Mapping
from typing import Any

import numpy as np

from dlamb import _core

# Copies of gstrickd constants (no runtime dependency on gstrickd). The tests in
# python/tests/test_gstrickd_contract.py compare them with the compiled gstrickd headers.

#: Largest number of photons in one chunk (gstrickd photon chunk limit).
CHUNK_ROWS = 1 << 16

#: Frames a source can be emitted in (gstrickd frame codes; G is not a source frame).
FRAMES = {"sky": 0, "telescope": 2, "camera": 3}

#: Columns of the emitted photons, in gstrickd column names.
COLUMNS = ("x", "y", "z", "cx", "cy", "cz", "t", "w", "wavelength")

# kind -> parameter names (in the order passed to the native core).
_SPECTRUM = {
    "fixed": ("wavelength_nm",),
    "uniform": ("low_nm", "high_nm"),
    "inverse_square": ("low_nm", "high_nm"),
    "table": ("wavelength_nm", "density"),
}
_PULSE = {
    "fixed": ("time_ns",),
    "tophat": ("width_ns", "delay_ns"),
    "gauss": ("sigma_ns", "mean_ns"),
    "exponential": ("decay_ns", "mean_ns"),
    "complex": ("flat_ns", "gauss_rms_ns", "decay_ns", "delay_ns"),
    "table": ("time_ns", "rate"),
}
_ANGULAR = {
    "isotropic": (),
    "cone": ("angle_deg",),
    "filled_cone": ("angle_deg",),
    "lambertian": (),
    "gauss": ("sigma_deg",),
    "rayleigh": (),
    "table": ("angle_deg", "value"),
}
_TARGET_KEYS = {"disc": ("centre_m", "normal", "radius_m"), "sphere": ("centre_m", "radius_m")}
_EMITTER_KEYS = ("position_m", "pointing", "angular")
_STAR_KEYS = ("towards_source", "distance_m", "disk_radius_m", "disk_z_m", "start_z_m")
_TOP_KEYS = ("frame", "photons", "emitter", "star", "spectrum", "pulse", "target")


def _check_keys(
    section: Mapping[str, Any], allowed: tuple[str, ...], label: str, *, has_kind: bool = False
) -> None:
    """Reject unknown keys; 'kind' is allowed only in profile and target sections."""
    extra = {"kind"} if has_kind else set()
    unknown = sorted(set(section) - set(allowed) - extra)
    if unknown:
        raise ValueError(f"unknown {label} key(s): {', '.join(unknown)}")


def _kind(section: Mapping[str, Any], table: Mapping[str, tuple[str, ...]], label: str) -> str:
    if "kind" not in section:
        raise ValueError(f"the {label} needs a 'kind' ({', '.join(table)})")
    kind = str(section["kind"])
    if kind not in table:
        raise ValueError(f"unknown {label} kind '{kind}' (one of {', '.join(table)})")
    return kind


# Optional parameters and their defaults (delays and mean times are zero).
_OPTIONAL = {"delay_ns": 0.0, "mean_ns": 0.0}


def _numbers(section: Mapping[str, Any], names: tuple[str, ...]) -> list[float]:
    values = []
    for name in names:
        if name not in section and name not in _OPTIONAL:
            raise ValueError(f"missing '{name}'")
        values.append(float(section.get(name, _OPTIONAL.get(name))))
    return values


def _profile(section, table, label, x_key, y_key):
    """Native (kind, params, x, y) of a spectrum, pulse or angular section."""
    kind = _kind(section, table, label)
    names = table[kind]
    _check_keys(section, names + ((x_key, y_key) if kind == "table" else ()), label, has_kind=True)
    if kind == "table":
        return kind, [], list(section.get(x_key, [])), list(section.get(y_key, []))
    return kind, _numbers(section, names), [], []


def _angular(section: Mapping[str, Any]):
    kind, params, x, y = _profile(section, _ANGULAR, "angular profile", "angle_deg", "value")
    if kind == "table":
        # The core takes u = 1 - cos(theta) (a density per solid angle).
        degrees = np.asarray(x, dtype=float)
        if (
            degrees.ndim != 1
            or not np.all(np.isfinite(degrees))
            or np.any((degrees < 0.0) | (degrees > 180.0))
        ):
            raise ValueError("angle_deg values must be finite and in [0, 180]")
        x = (1.0 - np.cos(np.radians(degrees))).tolist()
    return _core.make_angular(kind, params, x, y)


def _spectrum(section: Mapping[str, Any]):
    kind, params, x, y = _profile(section, _SPECTRUM, "spectrum", "wavelength_nm", "density")
    return _core.make_spectrum(kind, params, x, y)


def _pulse(section: Mapping[str, Any]):
    kind, params, x, y = _profile(section, _PULSE, "pulse", "time_ns", "rate")
    return _core.make_pulse(kind, params, x, y)


def build_source(spec: Mapping[str, Any]) -> _core.Source:
    """Native source of a source dictionary (validated, see the module docstring)."""
    _check_keys(spec, _TOP_KEYS, "source")
    if ("emitter" in spec) == ("star" in spec):
        raise ValueError("a source needs exactly one of 'emitter' or 'star'")
    frame = str(spec.get("frame", "sky"))
    if frame not in FRAMES:
        raise ValueError(f"unknown frame '{frame}' (one of {', '.join(FRAMES)})")
    if "photons" not in spec:
        raise ValueError("a source needs 'photons'")
    if "emitter" in spec:
        section = spec["emitter"]
        _check_keys(section, _EMITTER_KEYS, "emitter")
        for key in _EMITTER_KEYS:
            if key not in section:
                raise ValueError(f"the emitter needs '{key}'")
        geometry = _core.Emitter(
            [float(v) for v in section["position_m"]],
            [float(v) for v in section["pointing"]],
            _angular(section["angular"]),
        )
    else:
        section = spec["star"]
        _check_keys(section, _STAR_KEYS, "star")
        for key in ("towards_source", "disk_radius_m", "disk_z_m", "start_z_m"):
            if key not in section:
                raise ValueError(f"the star needs '{key}'")
        towards = np.asarray(section["towards_source"], dtype=float)
        if towards.shape != (3,) or not np.linalg.norm(towards) > 0.0:
            raise ValueError("towards_source needs three components, not all zero")
        geometry = _core.Star(
            (towards / np.linalg.norm(towards)).tolist(),
            float(section.get("distance_m", math.inf)),
            float(section["disk_radius_m"]),
            float(section["disk_z_m"]),
            float(section["start_z_m"]),
        )
    source = _core.Source(
        FRAMES[frame],
        geometry,
        _spectrum(spec.get("spectrum", {"kind": "fixed", "wavelength_nm": 400.0})),
        _pulse(spec.get("pulse", {"kind": "fixed", "time_ns": 0.0})),
        float(spec["photons"]),
    )
    if "target" in spec:
        target = spec["target"]
        if "emitter" not in spec:
            raise ValueError("only emitters take a target")
        kind = _kind(target, _TARGET_KEYS, "target")
        _check_keys(target, _TARGET_KEYS[kind], "target", has_kind=True)
        for key in _TARGET_KEYS[kind]:
            if key not in target:
                raise ValueError(f"the target needs '{key}'")
        if kind == "disc":
            source.set_disc_target(
                [float(v) for v in target["centre_m"]],
                [float(v) for v in target["normal"]],
                float(target["radius_m"]),
            )
        else:
            source.set_sphere_target(
                [float(v) for v in target["centre_m"]], float(target["radius_m"])
            )
    source.validate()
    return source


def expected_photons(spec: Mapping[str, Any]) -> float:
    """Mean number of generated photons of a source (candidates for a disc target)."""
    return build_source(spec).expected_count()


def photon_count(spec: Mapping[str, Any], *, seed: int, event: int = 0, telescope: int = 0) -> int:
    """Number of photons of a source in one (event, telescope), keyed on the seed."""
    return build_source(spec).photon_count(seed, event, telescope)


def emit_chunk(
    spec: Mapping[str, Any],
    *,
    seed: int,
    event: int = 0,
    telescope: int = 0,
    first_index: int = 0,
    rows: int = CHUNK_ROWS,
) -> dict[str, np.ndarray]:
    """Photons first_index .. first_index + rows - 1 (at most 2^16) as float32 columns."""
    if not 0 < rows <= CHUNK_ROWS:
        raise ValueError(f"rows must be in 1..{CHUNK_ROWS}")
    if first_index < 0:
        raise ValueError("first_index must be >= 0")
    return build_source(spec).emit(seed, event, telescope, first_index, rows)


def iter_chunks(
    spec: Mapping[str, Any], *, seed: int, event: int = 0, telescope: int = 0
) -> Iterator[dict[str, np.ndarray]]:
    """All photons of a source in chunks of at most 2^16 rows (bounded memory)."""
    source = build_source(spec)
    count = source.photon_count(seed, event, telescope)
    for first in range(0, count, CHUNK_ROWS):
        rows = min(CHUNK_ROWS, count - first)
        yield source.emit(seed, event, telescope, first, rows)


def emit(
    spec: Mapping[str, Any], *, seed: int, event: int = 0, telescope: int = 0
) -> dict[str, np.ndarray]:
    """All photons of a source as float32 columns (gstrickd names; see ``COLUMNS``).

    The arrays hold every photon in memory; use ``iter_chunks`` for large counts.
    """
    chunks = list(iter_chunks(spec, seed=seed, event=event, telescope=telescope))
    if not chunks:
        return {name: np.empty(0, dtype=np.float32) for name in COLUMNS}
    return {name: np.concatenate([chunk[name] for chunk in chunks]) for name in chunks[0]}
