"""Cherenkov light of a muon-like track on a telescope sphere, sky frame.

Port of the old ``obdeect sources muon-track`` program (``git show HEAD:`` of its header
``src/sources/cpp/include/obdeect/sources/programs/muon-track.hpp``). Physics as before:

- Track from the start height down to the impact point (x, y, 0), direction
  ``(-sin z cos a, -sin z sin a, -cos z)`` for zenith z and azimuth a.
- Cherenkov light with the Frank-Tamm yield of ``_cherenkov``. A constant index gives one line
  for the whole track. An index table or exponential gives lines of steps of 10 m, each with the
  index at its midpoint (the old program's stepping).
- Emission times: light reaches the impact point at t = 0.
- The telescope is a sphere; ``target: auto`` samples the photons that reach it.

Changes from the old program:

- Index as n - 1 only (the old program also accepted n). Index is required (no default).
- Wavelengths are required. Mass, charge and beta keep the old meaning; ``--primary`` defaults
  to muon (mu-), beta to 1.
- Core scattering (``--core-scatter-radius``, ``--core-scatter-count``): the old program
  rejected it. Here ``count`` core positions are drawn uniformly in a disc of the given radius
  around the impact point (``_common.choice_rng(--seed)``); each position is one event entry
  with ``--events`` repeats.
- Not ported: multiple scattering (rejected by the old program too), atmosphere tables
  (``-I``, ``-A``), and bunch size (no effect in the old program).

Run description: one telescope, one event entry per core position with ``repeat = --events``.
The event metadata holds primary id, energy [GeV] and core position.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from obdeect.sources.presets import _cherenkov, _common

NAME = "muon"
HELP = "Cherenkov light of a muon-like track on a telescope sphere (sky frame)"

#: Step length of index profiles [m] (numerical, as in the old program).
PROFILE_STEP_M = 10.0

#: Primary names for ``--primary``.
PRIMARY_NAMES = {
    "gamma": 1,
    "positron": 2,
    "electron": 3,
    "muon": 6,
    "mu-": 6,
    "mu+": 5,
    "proton": 14,
}

#: CORSIKA particle id -> (charge, mass [GeV]) for the particles of the old program's table.
_PARTICLES = {
    3: (-1, 0.000511),
    2: (1, 0.000511),
    6: (-1, 0.105658),
    5: (1, 0.105658),
    8: (1, 0.13957),
    11: (1, 0.493677),
    14: (1, 0.938272),
    13: (0, 0.939565),
    25: (0, 0.939565),
    9: (-1, 0.13957),
    12: (-1, 0.493677),
    10: (0, 0.497648),
    16: (0, 0.497648),
    15: (-1, 0.938272),
    66: (0, 1e-11),
    67: (0, 1e-11),
    68: (0, 1e-11),
    69: (0, 1e-11),
    133: (0, 1e-11),
    134: (0, 1e-11),
    130: (1, 1.777),
    131: (-1, 1.777),
    7: (0, 0.134),
    17: (0, 0.54751),
    50: (0, 0.78265),
    51: (0, 0.7755),
    52: (1, 0.7755),
    53: (-1, 0.7755),
}


def primary_id(text: str) -> int:
    """CORSIKA id from a number, or from a primary name (see PRIMARY_NAMES)."""
    if text[:1].isdigit():
        try:
            return int(text)
        except ValueError:
            raise ValueError(f"invalid primary id '{text}'") from None
    if text not in PRIMARY_NAMES:
        raise ValueError(f"unknown primary name: {text}")
    return PRIMARY_NAMES[text]


def mass_charge(pid: int) -> tuple[float, float]:
    """Mass [GeV] and charge of a CORSIKA particle id (as the old program)."""
    if pid in _PARTICLES:
        charge, mass = _PARTICLES[pid]
        return mass, float(charge)
    if pid <= 1:
        return 0.0, 0.0
    if 1 < pid < 100 or 110 < pid < 199:
        return 1.2, 0.0
    if pid == 402:
        return 0.931845 * 4.0, 2.0
    return (0.931494 - 0.000511) * (pid // 100), float(pid % 100)


def core_positions(
    x: float, y: float, radius_m: float, count: int, seed: int
) -> list[tuple[float, float]]:
    """Core positions: ``count`` uniform in a disc of ``radius_m`` around (x, y).

    With radius 0 only one position (x, y) is allowed.
    """
    if radius_m < 0.0 or not math.isfinite(radius_m):
        raise ValueError("core scatter radius must be finite and >= 0 m")
    if count < 1:
        raise ValueError("core scatter count must be >= 1")
    if radius_m == 0.0:
        if count != 1:
            raise ValueError("--core-scatter-count > 1 needs --core-scatter-radius > 0")
        return [(x, y)]
    rng = _common.choice_rng(seed)
    u = rng.random((count, 2))
    r = radius_m * np.sqrt(u[:, 0])
    phi = 2.0 * math.pi * u[:, 1]
    return [(x + float(ri * np.cos(p)), y + float(ri * np.sin(p))) for ri, p in zip(r, phi)]


def add_arguments(parser: argparse.ArgumentParser) -> None:
    _common.add_run_arguments(parser)
    _cherenkov.add_wavelength_arguments(parser)
    _cherenkov.add_index_arguments(parser)

    telescope = parser.add_argument_group("telescope")
    telescope.add_argument("--telescope-radius", type=float, required=True, metavar="M")
    telescope.add_argument(
        "--telescope-position",
        type=_common.vector,
        default="0,0,0",
        metavar="X,Y,Z",
        help="sphere centre [m] (default: origin)",
    )

    track = parser.add_argument_group("track")
    track.add_argument("--primary", default="muon", metavar="ID|NAME")
    track.add_argument("--mass", type=float, metavar="GEV", help="default: from the primary")
    track.add_argument("--charge", type=float, metavar="Q", help="default: from the primary")
    speed = track.add_mutually_exclusive_group()
    speed.add_argument("--momentum", type=float, metavar="GEV")
    speed.add_argument("--beta", type=float, metavar="B", help="default 1")
    track.add_argument("--impact-x", type=float, default=0.0, metavar="M")
    track.add_argument("--impact-y", type=float, default=0.0, metavar="M")
    track.add_argument("--start-height", type=float, required=True, metavar="M")
    track.add_argument("--zenith", type=float, default=0.0, metavar="DEG")
    track.add_argument("--azimuth", type=float, default=0.0, metavar="DEG")

    core = parser.add_argument_group("core scattering")
    core.add_argument("--core-scatter-radius", type=float, default=0.0, metavar="M")
    core.add_argument("--core-scatter-count", type=int, default=1, metavar="N")


def _check(args: argparse.Namespace) -> None:
    values = {
        "observation level": args.observation_level,
        "telescope radius": args.telescope_radius,
        "impact x": args.impact_x,
        "impact y": args.impact_y,
        "start height": args.start_height,
        "azimuth": args.azimuth,
    }
    for name, value in values.items():
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite, got {value}")
    if args.telescope_radius <= 0.0:
        raise ValueError("telescope radius must be > 0 m")
    if args.start_height <= 0.0:
        raise ValueError("start height must be > 0 m")
    if not (0.0 <= args.zenith < 90.0):
        raise ValueError(f"zenith must be in [0, 90) deg, got {args.zenith}")
    if args.events < 1:
        raise ValueError("events must be >= 1")
    if not (200.0 <= args.wavelength_min <= args.wavelength_max <= 2000.0):
        raise ValueError("wavelengths must satisfy 200 <= min <= max <= 2000 nm")
    if args.momentum is not None and not (args.momentum > 0.0 and math.isfinite(args.momentum)):
        raise ValueError("momentum must be > 0 GeV")
    if args.charge is not None and (args.charge == 0.0 or not math.isfinite(args.charge)):
        raise ValueError("charge must be finite and non-zero")
    if args.mass is not None and not (0.0005 <= args.mass <= 100.0):
        raise ValueError("mass must lie in [0.0005, 100] GeV")
    if args.beta is not None and not (0.0 < args.beta <= 1.0):
        raise ValueError("beta must lie in (0, 1]")


def build(args: argparse.Namespace) -> dict:
    """Run description of the muon track (see the module docstring)."""
    _check(args)
    pid = primary_id(args.primary)
    mass, charge = mass_charge(pid)
    if args.mass is not None:
        mass = args.mass
    if args.charge is not None:
        charge = args.charge
    if args.momentum is not None:
        beta = args.momentum / math.sqrt(args.momentum**2 + mass**2)
    else:
        beta = 1.0 if args.beta is None else args.beta
    energy_gev = mass / math.sqrt(1.0 - beta * beta) if beta < 1.0 and mass > 0.0 else 1e-9

    index = _cherenkov.index_from_args(args)
    step = None if args.n_minus_one is not None else PROFILE_STEP_M
    zenith = math.radians(args.zenith)
    azimuth = math.radians(args.azimuth)
    direction = [
        -math.sin(zenith) * math.cos(azimuth),
        -math.sin(zenith) * math.sin(azimuth),
        -math.cos(zenith),
    ]
    length_m = args.start_height / math.cos(zenith)

    events = []
    cores = core_positions(
        args.impact_x,
        args.impact_y,
        args.core_scatter_radius,
        args.core_scatter_count,
        args.seed,
    )
    for core_x, core_y in cores:
        emissions = _cherenkov.cherenkov_lines(
            end_m=[core_x, core_y, 0.0],
            direction=direction,
            length_m=length_m,
            step_m=step,
            index=index,
            beta=beta,
            wavelength_min_nm=args.wavelength_min,
            wavelength_max_nm=args.wavelength_max,
            charge=charge,
        )
        if not emissions:
            raise ValueError(
                "no Cherenkov light along the track: cos(theta_c) = 1/(n beta) > 1 everywhere"
            )
        events.append({
            "repeat": args.events,
            "metadata": {
                "primary_id": pid,
                "energy_gev": float(energy_gev),
                "core_x_m": core_x,
                "core_y_m": core_y,
            },
            "emissions": emissions,
        })

    return {
        "run": _common.run_section(args, frame="sky", preset=NAME),
        "telescopes": [
            {
                "id": 1,
                "position_m": list(args.telescope_position),
                "radius_m": args.telescope_radius,
            }
        ],
        "events": events,
    }
