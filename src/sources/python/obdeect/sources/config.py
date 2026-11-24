"""Run description: load from TOML or JSON, validate and fill defaults.

.. code-block:: toml

    [run]
    seed = 1
    observation_level_m = 2150.0

    [[telescopes]]
    radius_m = 12.5

    [[events]]
    [[events.emissions]]
    kind = "point"
    photons = 1e6
    position_m = [0.0, 0.0, 1.0]
    spectrum = { kind = "fixed", wavelength_nm = 400.0 }

``normalize`` returns a new dict with every default filled in (see
docs/sources/spec/run-file.md). Errors are ``ValueError`` whose message starts
with the key path, e.g. ``events[0].emissions[1].angular.kind: unknown kind 'cones' ...``.
Physics range checks and the errors of the core (plane waves, spectra, ...) are left to the C++
core; its messages start with the key path of the emission where it is known.
"""

import json
import math
import tomllib
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from obdeect.sources import _core

# Sentinels: a key that must be given, and an optional key without default (omitted).
REQUIRED = object()
ABSENT = object()

#: Frames of the photon stream (from the compiled stream contract).
_FRAMES = tuple(_core.stream_contract()["frames"])

#: Largest run seed: the seed is stored as a signed 64-bit integer in the stream.
MAX_SEED = 2**63 - 1


# Pairs of arrays that must have equal length (checked wherever both keys are present).
_PAIRS = (
    ("wavelength_nm", "density"),
    ("time_ns", "rate"),
    ("angle_deg", "value"),
    ("distance_m", "density"),
)


def load(path: str | Path) -> dict[str, Any]:
    """Read a run description from a .toml or .json file (not normalized)."""
    path = Path(path)
    if path.suffix not in (".toml", ".json"):
        raise ValueError(f"{path}: unsupported suffix '{path.suffix}' (use .toml or .json)")
    with path.open("rb") as stream:
        return tomllib.load(stream) if path.suffix == ".toml" else json.load(stream)


def normalize(config: Mapping[str, Any]) -> dict[str, Any]:
    """Validate a run description and return its normalized form (input is not changed)."""
    return _section(_TOP)(config, "")


def normalize_emission(spec: Mapping[str, Any]) -> dict[str, Any]:
    """Validate one emission (keys as in an event's emissions) and return its normalized form."""
    return _kinded(_EMISSION)(spec, "emission")


# --- checkers: each takes (value, path) and returns the normalized value -------------------


def _fail(path: str, message: str):
    raise ValueError(f"{path}: {message}" if path else message)


def _join(path: str, key: str) -> str:
    return f"{path}.{key}" if path else str(key)


def _type_name(value: Any) -> str:
    return type(value).__name__


def _range_text(lo, hi, strict=False) -> str:
    if hi is not None:
        return f"must be in {lo}..{hi}"
    return f"must be {'>' if strict else '>='} {lo}"


def _int(lo=None, hi=None):
    def check(value, path):
        if isinstance(value, bool) or not isinstance(value, int):
            _fail(path, f"expected an int, got {_type_name(value)}")
        if (lo is not None and value < lo) or (hi is not None and value > hi):
            _fail(path, _range_text(lo, hi))
        return value

    return check


def _num(lo=None, strict=False):
    def check(value, path):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            _fail(path, f"expected a number, got {_type_name(value)}")
        number = float(value)
        if not math.isfinite(number):
            _fail(path, "must be finite")
        if lo is not None and (number <= lo if strict else number < lo):
            _fail(path, _range_text(lo, None, strict))
        return number

    return check


_any_num = _num()


def _vec3(value, path):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        _fail(path, "expected a list of 3 numbers")
    return [_any_num(v, f"{path}[{i}]") for i, v in enumerate(value)]


def _column(value, path):
    """Array of a table (finite numbers, at least two)."""
    if not isinstance(value, (list, tuple)) or len(value) < 2:
        _fail(path, "expected a list of at least 2 numbers")
    return [_any_num(v, f"{path}[{i}]") for i, v in enumerate(value)]


def _bool(value, path):
    if not isinstance(value, bool):
        _fail(path, f"expected true or false, got {_type_name(value)}")
    return value


def _choice(options):
    def check(value, path):
        if value not in options:
            _fail(path, f"must be one of {', '.join(options)}")
        return value

    return check


def _metadata_value(item, path):
    if isinstance(item, str):
        return item
    if isinstance(item, list):
        return [_any_num(x, f"{path}[{i}]") for i, x in enumerate(item)]
    if isinstance(item, bool) or not isinstance(item, (int, float)):
        _fail(path, f"expected int, float, str or list of numbers, got {_type_name(item)}")
    return item if isinstance(item, int) else _any_num(item, path)


def _metadata(value, path):
    return {k: _metadata_value(v, _join(path, k)) for k, v in _mapping(value, path).items()}


def _mapping(value, path):
    if not isinstance(value, Mapping):
        _fail(path, f"expected a table, got {_type_name(value)}")
    return value


def _listof(check):
    """A non-empty list; each item is checked with its path, e.g. ``events[0]``."""

    def run(value, path):
        if not isinstance(value, list):
            _fail(path, f"expected a list of tables, got {_type_name(value)}")
        if not value:
            _fail(path, "needs at least one entry")
        return [check(item, f"{path}[{i}]") for i, item in enumerate(value)]

    return run


def _section(schema):
    """A table with the given fields."""
    return lambda value, path: _fields(_mapping(value, path), path, schema)


def _kinded(table):
    """A table with a 'kind' (one of the keys of ``table``) and the fields of that kind."""

    def check(value, path):
        section = _mapping(value, path)
        if "kind" not in section:
            _fail(_join(path, "kind"), "required")
        kind = section["kind"]
        if not isinstance(kind, str) or kind not in table:
            _fail(
                _join(path, "kind"),
                f"unknown kind {kind!r} (one of {', '.join(table)})",
            )
        body = {k: v for k, v in section.items() if k != "kind"}
        return {"kind": kind, **_fields(body, path, table[kind])}

    return check


def _fields(section: Mapping[str, Any], path: str, schema: dict) -> dict[str, Any]:
    """Check a table against ``schema``: key -> (checker, default); fills the defaults."""
    for key in section:
        if key not in schema:
            _fail(_join(path, key), "unknown key")
    out = {}
    for key, (check, default) in schema.items():
        child = _join(path, key)
        if key in section:
            out[key] = check(section[key], child)
        elif default is REQUIRED:
            _fail(child, "required")
        elif default is not ABSENT:
            out[key] = check(default, child)
    for x, y in _PAIRS:
        if isinstance(out.get(x), list) and isinstance(out.get(y), list):
            if len(out[x]) != len(out[y]):
                _fail(path, f"{x} and {y} differ in length")
    return out


def _target(value, path):
    if isinstance(value, str):  # shorthand: "none" or "auto"
        value = {"kind": value}
    return _kinded(_TARGET)(value, path)


def _telescopes(value, path):
    if isinstance(value, list):  # default id: position in the list, counted from 1
        value = [{"id": i, **v} if isinstance(v, Mapping) else v for i, v in enumerate(value, 1)]
    telescopes = _listof(_section(_TELESCOPE))(value, path)
    seen = {}
    for i, telescope in enumerate(telescopes):
        tid = telescope["id"]
        if tid in seen:
            _fail(f"{path}[{i}].id", f"duplicate id {tid} (also {path}[{seen[tid]}])")
        seen[tid] = i
    return telescopes


# --- schemas (key -> (checker, default)) ------------------------------------------------------

_SPECTRUM = {
    "fixed": {"wavelength_nm": (_any_num, REQUIRED)},
    "uniform": {"low_nm": (_any_num, REQUIRED), "high_nm": (_any_num, REQUIRED)},
    "inverse_square": {"low_nm": (_any_num, REQUIRED), "high_nm": (_any_num, REQUIRED)},
    "table": {"wavelength_nm": (_column, REQUIRED), "density": (_column, REQUIRED)},
}
_PULSE = {
    "fixed": {"time_ns": (_any_num, REQUIRED)},
    "tophat": {"width_ns": (_any_num, REQUIRED), "delay_ns": (_any_num, 0.0)},
    "gauss": {"sigma_ns": (_any_num, REQUIRED), "mean_ns": (_any_num, 0.0)},
    "exponential": {"decay_ns": (_any_num, REQUIRED), "mean_ns": (_any_num, 0.0)},
    "complex": {
        "flat_ns": (_any_num, REQUIRED),
        "gauss_rms_ns": (_any_num, REQUIRED),
        "decay_ns": (_any_num, REQUIRED),
        "delay_ns": (_any_num, 0.0),
    },
    "table": {"time_ns": (_column, REQUIRED), "rate": (_column, REQUIRED)},
}
_ANGULAR = {
    "isotropic": {},
    "cone": {"angle_deg": (_any_num, REQUIRED)},
    "filled_cone": {"angle_deg": (_any_num, REQUIRED)},
    "lambertian": {},
    "gauss": {"sigma_deg": (_num(0, strict=True), REQUIRED)},
    "rayleigh": {},
    "table": {"angle_deg": (_column, REQUIRED), "value": (_column, REQUIRED)},
}
_TARGET = {
    "none": {},
    "auto": {},
    "disc": {
        "centre_m": (_vec3, REQUIRED),
        "normal": (_vec3, REQUIRED),
        "radius_m": (_num(0, strict=True), REQUIRED),
    },
    "sphere": {"centre_m": (_vec3, REQUIRED), "radius_m": (_num(0, strict=True), REQUIRED)},
}
_LINE_PROFILE = {
    "distance_m": (_column, REQUIRED),
    "density": (_column, REQUIRED),
}

# Fields every emission has, and the geometry fields of each emission kind.
_COMMON = {
    "photons": (_num(0), REQUIRED),
    "bunch_size": (_num(1), 1.0),
    "time_offset_ns": (_any_num, 0.0),
    "to_plane_z0": (_bool, False),
    "spectrum": (_kinded(_SPECTRUM), REQUIRED),
    "pulse": (_kinded(_PULSE), {"kind": "fixed", "time_ns": 0.0}),
    "target": (_target, "none"),
}
_EMISSION = {
    "point": {
        "position_m": (_vec3, REQUIRED),
        "axis": (_vec3, [0.0, 0.0, 1.0]),
        "angular": (_kinded(_ANGULAR), {"kind": "isotropic"}),
    },
    "line": {
        "start_m": (_vec3, REQUIRED),
        "direction": (_vec3, REQUIRED),
        "length_m": (_any_num, REQUIRED),
        "beta": (_any_num, 1.0),
        "angular": (_kinded(_ANGULAR), {"kind": "cone", "angle_deg": 0.0}),
        "profile": (_section(_LINE_PROFILE), ABSENT),
    },
    "plane_wave": {
        "towards_source": (_vec3, REQUIRED),
        "distance_m": (_any_num, ABSENT),
        "disc_radius_m": (_num(0, strict=True), REQUIRED),
        "disc_z_m": (_any_num, 0.0),
        "start_z_m": (_any_num, ABSENT),
        "start_distance_m": (_num(0), ABSENT),
    },
    "sky_background": {
        "radius_m": (_num(0, strict=True), REQUIRED),
        "centre_m": (_vec3, [0.0, 0.0, 0.0]),
        "axis": (_vec3, [0.0, 0.0, -1.0]),
        "viewcone_deg": (_any_num, 90.0),
    },
}
_EMISSION = {kind: {**_COMMON, **geometry} for kind, geometry in _EMISSION.items()}

_RUN = {
    "number": (_int(1), 1),
    "seed": (_int(0, MAX_SEED), 1),
    "observation_level_m": (_any_num, REQUIRED),
    "atmosphere_id": (_int(), -1),
    "frame": (_choice(_FRAMES), "sky"),
    "metadata": (_metadata, {}),
}
_TELESCOPE = {
    "id": (_int(1), REQUIRED),
    "position_m": (_vec3, [0.0, 0.0, 0.0]),
    "radius_m": (_num(0, strict=True), REQUIRED),
}
_EVENT = {
    "repeat": (_int(0), 1),
    "metadata": (_metadata, {}),
    "emissions": (_listof(_kinded(_EMISSION)), REQUIRED),
}
_TOP = {
    "run": (_section(_RUN), REQUIRED),
    "telescopes": (_telescopes, REQUIRED),
    "events": (_listof(_section(_EVENT)), REQUIRED),
}
