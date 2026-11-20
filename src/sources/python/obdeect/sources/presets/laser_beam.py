"""laser-beam: primary photons of a laser beam inside a telescope field of view.

One line emission in the sky frame. It starts where the beam enters the field-of-view cone of
the telescope, points along the beam, and has the length of the path inside the cone. Its
photons are ``photons_per_m * length``. Time grows s/c along the beam (beta = 1). Only primary
photons are generated; scattering happens in obdeect propagation.

Geometry: the telescope is at the origin of the sky frame (x north, y west, z up). Beam and
telescope pointing are zenith and astronomical azimuth in degrees (0 = north, 90 = east,
clockwise seen from above). The beam zenith is the direction of propagation.
"""

from __future__ import annotations

import argparse

from obdeect.sources.presets import _common
from obdeect.sources.presets._beam import (
    cone_path,
    half_angle,
    non_negative,
    positive,
    unit_vector,
)

NAME = "laser-beam"
HELP = "primary photons of a laser beam inside a telescope field of view (one line emission)"


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the laser-beam preset (no device defaults: all geometry is required)."""
    _common.add_run_arguments(parser)
    _common.add_profile_arguments(parser, angular=None)
    group = parser.add_argument_group("laser")
    group.add_argument(
        "--laser-position",
        type=_common.vector,
        required=True,
        metavar="X,Y,Z",
        help="laser position in the sky frame [m]",
    )
    group.add_argument(
        "--laser-zenith",
        type=float,
        required=True,
        metavar="DEG",
        help="zenith of the beam propagation direction [deg]",
    )
    group.add_argument(
        "--laser-azimuth",
        type=float,
        required=True,
        metavar="DEG",
        help="astronomical azimuth of the beam: clockwise from north [deg]",
    )
    group.add_argument(
        "--photons-per-m",
        type=float,
        required=True,
        metavar="N",
        help="primary photons per metre of beam inside the field of view "
        "(sampling choice, not physics)",
    )
    group.add_argument(
        "--bunch-size", type=float, default=1.0, metavar="F", help="photons per bunch (>= 1)"
    )
    group = parser.add_argument_group("telescope")
    group.add_argument(
        "--telescope-zenith",
        type=float,
        required=True,
        metavar="DEG",
        help="zenith of the telescope pointing (cone axis) [deg]",
    )
    group.add_argument(
        "--telescope-azimuth",
        type=float,
        required=True,
        metavar="DEG",
        help="astronomical azimuth of the telescope pointing: clockwise from north [deg]",
    )
    group.add_argument(
        "--telescope-radius", type=float, required=True, metavar="M", help="telescope radius [m]"
    )
    group.add_argument(
        "--fov-deg",
        type=float,
        required=True,
        metavar="DEG",
        help="full opening angle of the field of view, in (0, 90) [deg]",
    )


def build(args: argparse.Namespace) -> dict:
    """Run description with one line emission per event."""
    half = half_angle(args.fov_deg)
    radius = positive(args.telescope_radius, "telescope radius")
    photons_per_m = non_negative(args.photons_per_m, "photons per metre")
    bunch = positive(args.bunch_size, "bunch size")
    if bunch < 1.0:
        raise ValueError(f"bunch size must be >= 1 (got {args.bunch_size})")

    beam = unit_vector(args.laser_zenith, args.laser_azimuth)
    axis = unit_vector(args.telescope_zenith, args.telescope_azimuth)
    start_at = tuple(args.laser_position)
    enter, length = cone_path(start_at, beam, axis, half)
    start = [start_at[i] + beam[i] * enter for i in range(3)]

    emission = {
        "kind": "line",
        "start_m": start,
        "direction": list(beam),
        "length_m": length,
        "beta": 1.0,
        "angular": {"kind": "cone", "angle_deg": 0.0},
        "photons": photons_per_m * length,
        "bunch_size": bunch,
        "spectrum": _common.profile(args.spectrum, "spectrum"),
        "pulse": _common.profile(args.pulse, "pulse"),
        "target": "none",
    }
    return {
        "run": _common.run_section(args, frame="sky"),
        "telescopes": [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": radius}],
        "events": [{"repeat": args.events, "emissions": [emission]}],
    }
