"""Shared builder of the pulsed-laser presets (``lidar``, ``laser-beam``).

Model: a pulsed laser at a position above the observation level points in any direction. Its
pulse is one point emission of primary photons, ``E lambda / (h c)`` photons with pulse energy
``E``, wavelength ``lambda`` and exact SI constants. The emission starts at the laser, or further
along the beam: at ``--start-height`` for a non-horizontal beam, at ``--start-distance`` for a
horizontal one. The scattered light is not generated here; it is computed by propagation (single
scattering, Rayleigh and Mie, beam extinction, transmission)::

    obdeect sources run lidar ... -o primaries.pst
    obdeect propagation run --transmission T --site-altitude H --seed N --mode weighted \\
        --scattering S --scattering-points 64 -i primaries.pst -o scattered.pst

``--mode weighted`` is needed because the primary weights exceed 2048.

Time: t = 0 is the laser firing. A later start shifts time: ``time_offset_ns`` is the path from
the laser to the start divided by c (vacuum).

Limits:

- A non-horizontal beam (``--start-height``): its tilt from the vertical must stay within
  ``--max-tilt`` (default 60 deg), measured from nadir for a downward beam (zenith > 90 deg) and
  from the zenith for an upward beam (zenith < 90 deg). ``--start-height`` must lie on the beam
  between the observation level and the laser.
- A horizontal beam (zenith 90 deg): no tilt limit. The beam stays at the laser height; the start
  is ``--start-distance`` metres along the beam from the laser (default 0: at the laser). The
  primary photons are valid; the scattering chain above is plane-parallel and only for beams within
  ``--max-tilt``, so it is not valid for horizontal beams. Horizontal photons (cz = 0) need the
  sphere representation (docs/stream/format.md).

The beam direction is zenith and astronomical azimuth (clockwise from north; sky frame x north, y
west). Telescope spheres are checked against the direct beam: a beam that hits one is an error
unless ``--allow-direct-beam`` is given.
"""

from __future__ import annotations

import argparse
import math

from obdeect.sources.presets import _common
from obdeect.sources.presets._beam import non_negative, positive, unit_vector

#: Planck constant [J s] and speed of light [m/s], exact SI values.
PLANCK_J_S = 6.62607015e-34
SPEED_OF_LIGHT_M_S = 299792458.0
SPEED_OF_LIGHT_M_NS = SPEED_OF_LIGHT_M_S * 1e-9

Vector = tuple[float, float, float]


def photons_per_pulse(energy_j: float, wavelength_nm: float) -> float:
    """Photons in a pulse of energy ``E`` [J] at wavelength ``lambda`` [nm]: E lambda / (h c)."""
    return energy_j * (wavelength_nm * 1e-9) / (PLANCK_J_S * SPEED_OF_LIGHT_M_S)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Beam, pulse and telescope options shared by the pulsed-laser presets."""
    group = parser.add_argument_group("beam")
    group.add_argument(
        "--beam-zenith",
        type=float,
        required=True,
        metavar="DEG",
        help="zenith of the beam direction: > 90 points down, < 90 up, 90 is horizontal [deg]",
    )
    group.add_argument(
        "--beam-azimuth",
        type=float,
        required=True,
        metavar="DEG",
        help="astronomical azimuth of the beam: clockwise from north [deg]",
    )
    group.add_argument(
        "--max-tilt",
        type=float,
        default=60.0,
        metavar="DEG",
        help="largest tilt of a non-horizontal beam from the vertical in its own hemisphere (from "
        "nadir for a downward beam, from the zenith for an upward one), in [0, 90) (default 60); "
        "not used for a horizontal beam [deg]",
    )
    group.add_argument(
        "--divergence",
        type=float,
        required=True,
        metavar="DEG",
        help="full beam divergence angle, in [0, 180) (0: parallel beam) [deg]",
    )
    group.add_argument(
        "--start-height",
        type=float,
        metavar="M",
        help="non-horizontal beams: start the emission at this height on the beam, e.g. 100000 "
        "for a satellite (default: at the laser) [m]",
    )
    group.add_argument(
        "--start-distance",
        type=float,
        metavar="M",
        help="horizontal beams: start the emission this distance along the beam from the laser "
        "(default 0: at the laser) [m]",
    )
    group = parser.add_argument_group("pulse")
    group.add_argument(
        "--pulse-energy", type=float, required=True, metavar="J", help="pulse energy [J]"
    )
    group.add_argument(
        "--wavelength-nm",
        type=float,
        required=True,
        metavar="NM",
        help="wavelength [nm]; the spectrum is fixed at this wavelength",
    )
    group.add_argument(
        "--pulse",
        default="fixed:time_ns=0",
        metavar="KIND[:K=V,...]",
        help="time profile of the pulse (default fixed:time_ns=0)",
    )
    group.add_argument(
        "--rows",
        type=float,
        default=1e4,
        metavar="N",
        help="sampled photon rows per pulse; bunch size = photons / rows (default 1e4). "
        "A sampling choice: the weights carry the physical photons",
    )
    group = parser.add_argument_group("telescope")
    where = group.add_mutually_exclusive_group(required=True)
    where.add_argument(
        "--telescope-radius", type=float, metavar="M", help="one telescope at the origin [m]"
    )
    where.add_argument(
        "--telescope",
        type=_telescope,
        action="append",
        metavar="X,Y,Z,R",
        help="telescope position and radius in the sky frame [m]; repeatable (ids 1, 2, ...)",
    )
    group.add_argument(
        "--allow-direct-beam",
        action="store_true",
        help="keep the run when the direct laser beam hits a telescope sphere",
    )


def _telescope(text: str) -> dict:
    items = text.split(",")
    try:
        x, y, z, r = (float(v) for v in items)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{text}' needs four numbers X,Y,Z,R") from None
    return {"position_m": [x, y, z], "radius_m": r}


def direct_hit(
    start: Vector, beam: Vector, centre: Vector, radius: float, tan_half: float
) -> tuple[float, float, float] | None:
    """Direct-beam check of one telescope sphere.

    The beam is the cone of half angle ``atan(tan_half)`` about the ray ``start + t beam``, t >= 0.
    Returns ``(along, distance, beam_radius)`` at the closest approach to the centre if the cone
    reaches the sphere (``distance < radius + beam_radius``), else None. The closest point must lie
    at or below the start (t >= 0).
    """
    rel = [centre[i] - start[i] for i in range(3)]
    along = sum(rel[i] * beam[i] for i in range(3))
    if along < 0.0:
        return None
    perp = [rel[i] - along * beam[i] for i in range(3)]
    distance = math.sqrt(sum(p * p for p in perp))
    beam_radius = along * tan_half
    if distance < radius + beam_radius:
        return along, distance, beam_radius
    return None


def _start_on_beam(position, beam, height) -> tuple[tuple[float, float, float], float]:
    """Start of a non-horizontal beam at ``height`` and the time offset of the path from the laser.

    The start lies on the beam beyond the laser's observation side: below the laser for a downward
    beam, above it for an upward one; never below the observation level. No height: the laser.
    """
    if height is None:
        return tuple(position), 0.0
    along = (height - position[2]) * beam[2]
    if not (math.isfinite(height) and height >= 0.0 and along > 0.0):
        raise ValueError(
            f"--start-height {height} m must lie on the beam between the observation level "
            f"and the laser at {position[2]:g} m (below the laser for a downward beam, above "
            "it for an upward one)"
        )
    path = (height - position[2]) / beam[2]  # > 0
    start = tuple(position[i] + beam[i] * path for i in range(2)) + (height,)
    return start, path / SPEED_OF_LIGHT_M_NS


def build(args: argparse.Namespace, position) -> dict:
    """Run description of one pulse: one point emission at the laser or further along the beam.

    ``position`` is the laser position in the sky frame, z above the observation level.
    """
    position = tuple(float(v) for v in position)
    energy = positive(args.pulse_energy, "pulse energy")
    wavelength = positive(args.wavelength_nm, "wavelength")
    rows = positive(args.rows, "rows")
    divergence = non_negative(args.divergence, "divergence")
    if divergence >= 180.0:
        raise ValueError(f"divergence {divergence} deg must be < 180")
    max_tilt = non_negative(args.max_tilt, "max tilt")
    if max_tilt >= 90.0:
        raise ValueError(f"--max-tilt {max_tilt} deg must be < 90")
    if position[2] <= 0.0:
        raise ValueError(f"laser height {position[2]} m must be above the observation level")

    beam = unit_vector(args.beam_zenith, args.beam_azimuth)
    horizontal = abs(beam[2]) <= 1e-12
    start, time_offset_ns = position, 0.0
    if horizontal:
        if args.start_height is not None:
            raise ValueError(
                "--start-height is for non-horizontal beams: a horizontal beam stays at the laser "
                "height, so give --start-distance (metres along the beam from the laser)"
            )
        distance = (
            0.0
            if args.start_distance is None
            else non_negative(args.start_distance, "start distance")
        )
        start = tuple(position[i] + beam[i] * distance for i in range(3))
        time_offset_ns = distance / SPEED_OF_LIGHT_M_NS
    else:
        if args.start_distance is not None:
            raise ValueError(
                "--start-distance is for horizontal beams; a non-horizontal beam starts at "
                "--start-height"
            )
        tilt = math.degrees(math.acos(min(1.0, abs(beam[2]))))
        if tilt > max_tilt + 1e-9:  # the limit is inclusive; 1e-9 deg absorbs rounding
            side = "nadir" if beam[2] < 0.0 else "the zenith"
            raise ValueError(
                f"beam tilted {tilt:.3g} deg from {side} exceeds --max-tilt {max_tilt:g} deg: "
                "the scattering model is plane-parallel, so only near-vertical beams are valid"
            )
        start, time_offset_ns = _start_on_beam(position, beam, args.start_height)
    if args.telescope is not None:
        telescopes = args.telescope
    else:
        telescopes = [{"position_m": [0.0, 0.0, 0.0], "radius_m": args.telescope_radius}]
    for number, telescope in enumerate(telescopes, start=1):
        positive(telescope["radius_m"], f"telescope {number} radius")

    tan_half = math.tan(math.radians(0.5 * divergence))
    if not args.allow_direct_beam:
        for number, telescope in enumerate(telescopes, start=1):
            hit = direct_hit(start, beam, telescope["position_m"], telescope["radius_m"], tan_half)
            if hit is not None:
                along, distance, beam_radius = hit
                raise ValueError(
                    f"the direct laser beam hits telescope {number}: its closest approach is "
                    f"{distance:.3g} m from the centre, {along:.3g} m along the beam, where the "
                    f"beam radius is {beam_radius:.3g} m (telescope radius "
                    f"{telescope['radius_m']:g} m). Move the telescope or the beam, or pass "
                    "--allow-direct-beam"
                )

    per_pulse = photons_per_pulse(energy, wavelength)
    if divergence == 0.0:
        angular = {"kind": "cone", "angle_deg": 0.0}
    else:
        angular = {"kind": "filled_cone", "angle_deg": 0.5 * divergence}
    emission = {
        "kind": "point",
        "position_m": list(start),
        "axis": list(beam),
        "angular": angular,
        "photons": per_pulse,
        "bunch_size": max(1.0, per_pulse / rows),
        "time_offset_ns": time_offset_ns,
        "to_plane_z0": False,
        "spectrum": {"kind": "fixed", "wavelength_nm": wavelength},
        "pulse": _common.profile(args.pulse, "pulse"),
        "target": "none",
    }
    metadata = {
        "pulse_energy_j": energy,
        "wavelength_nm": wavelength,
        "divergence_deg": divergence,
        "beam_zenith_deg": args.beam_zenith,
        "beam_azimuth_deg": args.beam_azimuth,
        "laser_x_m": position[0],
        "laser_y_m": position[1],
        "laser_z_m": position[2],
        "start_height_m": start[2],
        "start_distance_m": math.dist(start, position),
    }
    return {
        "run": _common.run_section(args, frame="sky", **metadata),
        "telescopes": [
            {"id": i, "position_m": list(t["position_m"]), "radius_m": t["radius_m"]}
            for i, t in enumerate(telescopes, start=1)
        ],
        "events": [{"repeat": args.events, "emissions": [emission]}],
    }
