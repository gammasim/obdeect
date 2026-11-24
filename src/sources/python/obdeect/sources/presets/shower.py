"""Simple shower model: Cherenkov and fluorescence light of a straight air shower, sky frame.

Model (the spirit of the ``showermodel`` package, without its shower generation):

- One straight shower axis from the start height down to the core (core_x, core_y, 0). The
  light comes from the direction of ``--zenith`` and ``--azimuth`` (astronomical, sky frame; see
  presets.md), so the axis runs along their reverse.
- Longitudinal profile, Gaisser-Hillas in slant depth X [g/cm^2]:
  ``N(X) = Nmax ((X - X0)/(Xmax - X0))^((Xmax - X0)/lambda) exp((Xmax - X)/lambda)`` for X > X0,
  zero otherwise. X0 is the first-interaction depth, so a start height above it is allowed.
  No Greisen parameterisation: Nmax, Xmax, X0 and lambda are inputs.
- Atmosphere: vertical depth ``X_v(h) = X_vert0 exp(-h / H)`` with h above the observation level.
  The slant depth along the axis is ``X = X_v(h) / cos(zenith)``.
- The axis is split into equal steps of at most ``--step``. Per step, with N the particle number
  at the midpoint depth:
  - Cherenkov: one line. Electrons with beta = 1, yield from ``_cherenkov`` (Frank-Tamm),
    photons = yield x N x step, cone at the Cherenkov angle from n at the step midpoint,
    inverse-square spectrum between the two wavelength options.
  - Fluorescence: one isotropic line, photons = ``--fluorescence-yield`` x N x step, spectrum
    fixed or table (``--fluorescence-spectrum``).
- Emission times: the light front moves at c towards the core, t = 0 at the core.

Limits:

- No angular spread of the electrons: all Cherenkov light of a step is on the cone of the axis.
- No lateral distribution: all light is emitted on the axis.
- No multiple scattering.
- The Cherenkov angle uses n at the step midpoint (n is constant within a step).
- Electrons have beta = 1 and no energy dependence.
- Fluorescence yield per particle and metre is an input, not a model.

Run description: one telescope, one event entry with ``repeat = --events``. The event metadata
holds the core position and Xmax.
"""

from __future__ import annotations

import argparse
import math
from collections.abc import Callable

from obdeect.sources.presets import _cherenkov, _common
from obdeect.sources.presets._beam import travel_from_source

NAME = "shower"
HELP = "Cherenkov and fluorescence light of a straight air shower (simple model, sky frame)"


def gaisser_hillas(
    depth_g_cm2: float, nmax: float, xmax_g_cm2: float, x0_g_cm2: float, lambda_g_cm2: float
) -> float:
    """Gaisser-Hillas particle number at slant depth X; zero for X <= X0."""
    if depth_g_cm2 <= x0_g_cm2:
        return 0.0
    ratio = (depth_g_cm2 - x0_g_cm2) / (xmax_g_cm2 - x0_g_cm2)
    exponent = (xmax_g_cm2 - x0_g_cm2) / lambda_g_cm2
    return nmax * ratio**exponent * math.exp((xmax_g_cm2 - depth_g_cm2) / lambda_g_cm2)


def slant_depth(altitude_m: float, depth0_g_cm2: float, scale_height_m: float, zenith_deg: float):
    """Slant depth at altitude h: X_vert0 exp(-h / H) / cos(zenith)."""
    return (
        depth0_g_cm2 * math.exp(-altitude_m / scale_height_m) / math.cos(math.radians(zenith_deg))
    )


def _positive(value: float, name: str) -> None:
    if not (math.isfinite(value) and value > 0.0):
        raise ValueError(f"{name} must be finite and > 0, got {value}")


def add_arguments(parser: argparse.ArgumentParser) -> None:
    _cherenkov.add_shared_arguments(parser)

    axis = parser.add_argument_group("shower axis")
    axis.add_argument("--start-height", type=float, required=True, metavar="M")
    axis.add_argument(
        "--zenith",
        type=float,
        required=True,
        metavar="DEG",
        help="zenith of the direction the light comes from (the shower runs the other way)",
    )
    axis.add_argument(
        "--azimuth",
        type=float,
        required=True,
        metavar="DEG",
        help="astronomical azimuth of the light source (clockwise from north)",
    )
    axis.add_argument("--core-x", type=float, default=0.0, metavar="M")
    axis.add_argument("--core-y", type=float, default=0.0, metavar="M")
    axis.add_argument("--step", type=float, required=True, metavar="M", help="step length")

    profile = parser.add_argument_group("longitudinal profile (Gaisser-Hillas)")
    profile.add_argument("--nmax", type=float, required=True, metavar="N")
    profile.add_argument("--xmax", type=float, required=True, metavar="G_CM2")
    profile.add_argument("--x0", type=float, required=True, metavar="G_CM2")
    profile.add_argument("--gh-lambda", type=float, required=True, metavar="G_CM2")

    atmosphere = parser.add_argument_group("atmosphere")
    atmosphere.add_argument(
        "--atmosphere-depth0",
        type=float,
        required=True,
        metavar="G_CM2",
        help="vertical depth at the observation level",
    )
    atmosphere.add_argument("--scale-height", type=float, required=True, metavar="M")

    fluorescence = parser.add_argument_group("fluorescence")
    fluorescence.add_argument(
        "--fluorescence-yield",
        type=float,
        required=True,
        metavar="PHOTONS_PER_M",
        help="photons per metre and particle, into 4 pi",
    )
    fluorescence.add_argument(
        "--fluorescence-spectrum",
        required=True,
        metavar="fixed:K=V | table:FILE",
        help="fixed:wavelength_nm=337 or table:FILE (wavelength_nm density)",
    )


def _check(args: argparse.Namespace) -> None:
    for name, value in {
        "observation level": args.observation_level,
        "telescope radius": args.telescope_radius,
        "start height": args.start_height,
        "core x": args.core_x,
        "core y": args.core_y,
        "azimuth": args.azimuth,
    }.items():
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite, got {value}")
    _positive(args.telescope_radius, "telescope radius")
    _positive(args.start_height, "start height")
    _positive(args.step, "step length")
    _positive(args.nmax, "Nmax")
    _positive(args.gh_lambda, "lambda")
    _positive(args.atmosphere_depth0, "atmosphere depth")
    _positive(args.scale_height, "scale height")
    _positive(args.fluorescence_yield, "fluorescence yield")
    if not (0.0 <= args.zenith < 90.0):
        raise ValueError(f"zenith must be in [0, 90) deg, got {args.zenith}")
    if not (args.x0 < args.xmax and math.isfinite(args.x0) and math.isfinite(args.xmax)):
        raise ValueError("Gaisser-Hillas needs X0 < Xmax")
    if args.events < 1:
        raise ValueError("events must be >= 1")


def _fluorescence_spectrum(text: str) -> dict:
    spectrum = _common.profile(text, "spectrum")
    if spectrum["kind"] not in ("fixed", "table"):
        raise ValueError("fluorescence spectrum must be fixed or table")
    return spectrum


def fluorescence_lines(
    *,
    end_m,
    direction,
    length_m: float,
    step_m: float,
    yield_per_m: float,
    spectrum: dict,
    particles: Callable[[float], float],
) -> list[dict]:
    """Isotropic line emissions of fluorescence, one per step, with beta = 1."""
    lines = []
    for seg in _cherenkov.track_segments(end_m, direction, length_m, step_m):
        photons = yield_per_m * seg.length_m * particles(seg.altitude_m)
        if photons <= 0.0:
            continue
        lines.append({
            "kind": "line",
            "start_m": list(seg.start_m),
            "direction": [float(v) for v in direction],
            "length_m": seg.length_m,
            "beta": 1.0,
            "photons": photons,
            "time_offset_ns": _cherenkov.time_offset_ns(seg.distance_m, 1.0),
            "angular": {"kind": "isotropic"},
            "spectrum": dict(spectrum),
            "target": "auto",
        })
    return lines


def build(args: argparse.Namespace) -> dict:
    """Run description of the simple shower (see the module docstring)."""
    _check(args)
    index = _cherenkov.index_from_args(args)
    direction = list(travel_from_source(args.zenith, args.azimuth))
    length_m = args.start_height / math.cos(math.radians(args.zenith))
    end = [args.core_x, args.core_y, 0.0]

    def particles(altitude_m: float) -> float:
        depth = slant_depth(altitude_m, args.atmosphere_depth0, args.scale_height, args.zenith)
        return gaisser_hillas(depth, args.nmax, args.xmax, args.x0, args.gh_lambda)

    cherenkov = _cherenkov.cherenkov_lines(
        end_m=end,
        direction=direction,
        length_m=length_m,
        step_m=args.step,
        index=index,
        beta=1.0,
        wavelength_min_nm=args.wavelength_min,
        wavelength_max_nm=args.wavelength_max,
        charge=1.0,
        particles=particles,
    )
    fluorescence = fluorescence_lines(
        end_m=end,
        direction=direction,
        length_m=length_m,
        step_m=args.step,
        yield_per_m=args.fluorescence_yield,
        spectrum=_fluorescence_spectrum(args.fluorescence_spectrum),
        particles=particles,
    )
    emissions = cherenkov + fluorescence
    if not emissions:
        raise ValueError("no light: the profile is zero along the whole axis")

    return {
        "run": _common.run_section(args, frame="sky", preset=NAME),
        "telescopes": _common.single_telescope(args.telescope_radius, args.telescope_position),
        "events": [
            {
                "repeat": args.events,
                "metadata": {
                    "core_x_m": args.core_x,
                    "core_y_m": args.core_y,
                    "xmax_g_cm2": args.xmax,
                },
                "emissions": emissions,
            }
        ],
    }
