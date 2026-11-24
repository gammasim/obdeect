"""Cherenkov light of a charged particle on a straight track, shared by the muon and shower presets.

Conventions: lengths in m; altitudes ``z`` in m above the observation level (z up); wavelengths
in nm; angles in deg; times in ns. A track ends at a point ``end_m`` (the impact point or the
shower core) where the arrival time is t = 0. Light emitted at distance ``d`` before that point,
along the track, gets ``time_offset_ns = -d / (beta c)``.

Refractive index n - 1 against altitude, one of:

- a constant ``eta``;
- an exponential ``eta0 exp(-h / scale_height)``;
- a table file (lines ``altitude_m n-1``), linear in altitude, constant beyond the ends.

Frank-Tamm yield, mean photons per metre into 4 pi between wavelengths l1 < l2:

    dN/dx = 2 pi alpha z^2 sin^2(theta) (1/l1 - 1/l2) 1e9,   cos(theta) = 1 / (n beta)

A track is split into equal steps of at most ``step_m``. Each step is one line emission with a
cone at the Cherenkov angle from n at its midpoint, photons = yield x step x particles(z_mid),
and an inverse-square spectrum between l1 and l2. Steps where cos(theta) > 1 emit nothing.
"""

from __future__ import annotations

import argparse
import math
from collections.abc import Callable
from typing import NamedTuple

import numpy as np
from obdeect.sources.presets import _common
from obdeect.sources.presets._common import SPEED_OF_LIGHT_M_PER_NS, read_table

#: Fine-structure constant.
ALPHA = 0.0072973526

#: Upper limit of n - 1 accepted for any profile (air is about 3e-4).
MAX_ETA = 0.01

IndexFunction = Callable[[float], float]


class Segment(NamedTuple):
    """One step of a track: upper end, length, distance from the end point, midpoint altitude."""

    start_m: tuple[float, float, float]
    length_m: float
    distance_m: float
    altitude_m: float


def _check_eta(eta: float, name: str) -> None:
    if not (math.isfinite(eta) and 0.0 < eta < MAX_ETA):
        raise ValueError(f"{name} must be in (0, {MAX_ETA}), got {eta}")


def constant_index(eta: float) -> IndexFunction:
    """Refractive index n - 1 = eta at all altitudes."""
    _check_eta(eta, "n - 1")
    return lambda h: eta


def exponential_index(eta0: float, scale_height_m: float) -> IndexFunction:
    """Refractive index n - 1 = eta0 exp(-h / scale_height_m)."""
    _check_eta(eta0, "n - 1 at h = 0")
    if not (math.isfinite(scale_height_m) and scale_height_m > 0.0):
        raise ValueError(f"scale height must be finite and > 0 m, got {scale_height_m}")
    return lambda h: eta0 * math.exp(-h / scale_height_m)


def table_index(altitude_m, eta) -> IndexFunction:
    """Refractive index linear in altitude, constant beyond the first and last node."""
    x = np.asarray(altitude_m, dtype=float)
    y = np.asarray(eta, dtype=float)
    if x.ndim != 1 or x.size == 0 or x.size != y.size:
        raise ValueError("index table needs equal numbers of altitudes and n - 1 values")
    if not np.all(np.isfinite(x)) or np.any(np.diff(x) <= 0.0):
        raise ValueError("index table altitudes must be finite and strictly increasing")
    if not np.all(np.isfinite(y)) or np.any((y <= 0.0) | (y >= MAX_ETA)):
        raise ValueError(f"index table n - 1 values must be in (0, {MAX_ETA})")
    return lambda h: float(np.interp(h, x, y))


def index_from_file(path: str) -> IndexFunction:
    """Refractive index from a two-column table file (altitude_m, n - 1)."""
    altitude, eta = read_table(path)
    return table_index(altitude, eta)


def cherenkov_cos(eta: float, beta: float) -> float:
    """cos(theta_c) = 1 / (n beta). Values above 1 mean no Cherenkov light."""
    if not (math.isfinite(beta) and 0.0 < beta <= 1.0):
        raise ValueError(f"beta must be in (0, 1], got {beta}")
    return 1.0 / ((1.0 + eta) * beta)


def _check_wavelengths(wavelength_min_nm: float, wavelength_max_nm: float) -> None:
    if not (0.0 < wavelength_min_nm < wavelength_max_nm and math.isfinite(wavelength_max_nm)):
        raise ValueError(
            f"wavelengths must satisfy 0 < min < max, got {wavelength_min_nm}, {wavelength_max_nm}"
        )


def frank_tamm_yield(
    cos_theta: float, wavelength_min_nm: float, wavelength_max_nm: float, charge: float = 1.0
) -> float:
    """Mean Cherenkov photons per metre (into 4 pi) of charge ``charge`` at cos(theta)."""
    if not (0.0 < cos_theta <= 1.0):
        raise ValueError(f"cos(theta) must be in (0, 1], got {cos_theta}")
    _check_wavelengths(wavelength_min_nm, wavelength_max_nm)
    sin2 = 1.0 - cos_theta * cos_theta
    band = 1.0 / wavelength_min_nm - 1.0 / wavelength_max_nm
    return 2.0 * math.pi * ALPHA * charge * charge * sin2 * band * 1e9


def time_offset_ns(distance_m: float, beta: float) -> float:
    """Emission time of light emitted ``distance_m`` before the end point (arrival at t = 0)."""
    if not (math.isfinite(beta) and 0.0 < beta <= 1.0):
        raise ValueError(f"beta must be in (0, 1], got {beta}")
    return -distance_m / (beta * SPEED_OF_LIGHT_M_PER_NS)


def track_segments(end_m, direction, length_m: float, step_m: float | None = None) -> list[Segment]:
    """Equal steps of a straight track of ``length_m`` ending at ``end_m``, lowest step first.

    ``direction`` is the unit direction of motion (downward for a track from above). With
    ``step_m=None`` the track is one step. The number of steps is ceil(length / step).
    """
    end = np.asarray(end_m, dtype=float)
    unit = np.asarray(direction, dtype=float)
    if end.shape != (3,) or not np.all(np.isfinite(end)):
        raise ValueError("track end must be three finite numbers")
    if unit.shape != (3,) or not math.isclose(float(np.linalg.norm(unit)), 1.0, abs_tol=1e-9):
        raise ValueError("track direction must be a unit vector")
    if not (math.isfinite(length_m) and length_m > 0.0):
        raise ValueError(f"track length must be finite and > 0 m, got {length_m}")
    if step_m is None:
        count = 1
    else:
        if not (math.isfinite(step_m) and step_m > 0.0):
            raise ValueError(f"step length must be finite and > 0 m, got {step_m}")
        count = max(1, math.ceil(length_m / step_m))
    step = length_m / count
    upward = -unit
    segments = []
    for k in range(count):
        upper = (k + 1) * step
        start = end + upward * upper
        midpoint = end + upward * ((k + 0.5) * step)
        segments.append(
            Segment(
                start_m=(float(start[0]), float(start[1]), float(start[2])),
                length_m=step,
                distance_m=upper,
                altitude_m=float(midpoint[2]),
            )
        )
    return segments


def cherenkov_lines(
    *,
    end_m,
    direction,
    length_m: float,
    step_m: float | None,
    index: IndexFunction,
    beta: float,
    wavelength_min_nm: float,
    wavelength_max_nm: float,
    charge: float = 1.0,
    particles: Callable[[float], float] | None = None,
) -> list[dict]:
    """Cherenkov line emissions of a straight track (run description emissions).

    ``particles(z)`` gives the number of particles per step at midpoint altitude z (default 1).
    Steps without Cherenkov light or without photons are left out.
    """
    _check_wavelengths(wavelength_min_nm, wavelength_max_nm)
    lines = []
    for seg in track_segments(end_m, direction, length_m, step_m):
        cos_theta = cherenkov_cos(index(seg.altitude_m), beta)
        if cos_theta > 1.0:
            continue
        number = 1.0 if particles is None else particles(seg.altitude_m)
        photons = (
            frank_tamm_yield(cos_theta, wavelength_min_nm, wavelength_max_nm, charge)
            * seg.length_m
            * number
        )
        if photons <= 0.0:
            continue
        lines.append({
            "kind": "line",
            "start_m": list(seg.start_m),
            "direction": [float(v) for v in direction],
            "length_m": seg.length_m,
            "beta": float(beta),
            "photons": photons,
            "time_offset_ns": time_offset_ns(seg.distance_m, beta),
            "angular": {"kind": "cone", "angle_deg": math.degrees(math.acos(cos_theta))},
            "spectrum": {
                "kind": "inverse_square",
                "low_nm": float(wavelength_min_nm),
                "high_nm": float(wavelength_max_nm),
            },
            "target": "auto",
        })
    return lines


def add_index_arguments(parser: argparse.ArgumentParser) -> None:
    """Exactly one refractive index option (constant, exponential or table)."""
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--n-minus-one", type=float, metavar="ETA", help="constant n - 1")
    group.add_argument(
        "--index-exponential",
        nargs=2,
        type=float,
        metavar=("ETA0", "SCALE_HEIGHT_M"),
        help="n - 1 = ETA0 exp(-h / SCALE_HEIGHT_M)",
    )
    group.add_argument(
        "--index-profile", metavar="FILE", help="table: lines 'altitude_m n-1' (# comments)"
    )


def index_from_args(args: argparse.Namespace) -> IndexFunction:
    """Refractive index function from the options of :func:`add_index_arguments`."""
    if args.n_minus_one is not None:
        return constant_index(args.n_minus_one)
    if args.index_exponential is not None:
        eta0, scale_height_m = args.index_exponential
        return exponential_index(eta0, scale_height_m)
    return index_from_file(args.index_profile)


def add_shared_arguments(parser: argparse.ArgumentParser) -> None:
    """Options shared by the Cherenkov presets (muon, shower): run, band, index and telescope."""
    _common.add_run_arguments(parser)
    add_wavelength_arguments(parser)
    add_index_arguments(parser)
    telescope = parser.add_argument_group("telescope")
    telescope.add_argument("--telescope-radius", type=float, required=True, metavar="M")
    telescope.add_argument(
        "--telescope-position",
        type=_common.vector,
        default="0,0,0",
        metavar="X,Y,Z",
        help="sphere centre [m] (default: origin)",
    )


def add_wavelength_arguments(parser: argparse.ArgumentParser) -> None:
    """Wavelength band of the Cherenkov spectrum (inverse square between the two limits)."""
    group = parser.add_argument_group("Cherenkov band")
    group.add_argument("--wavelength-min", type=float, required=True, metavar="NM")
    group.add_argument("--wavelength-max", type=float, required=True, metavar="NM")
