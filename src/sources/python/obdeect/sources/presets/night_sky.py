"""night-sky: night-sky background photons on the sphere of one telescope (sky frame).

One ``sky_background`` emission per event. It is centred on the telescope with the telescope
radius, and its axis is minus the telescope pointing (the telescope looks up at the sky, so
the axis points down). ``viewcone-deg`` is the half opening angle of the viewed cone. Photons
per event are ``--photons``. The spectrum is required, for example ``table:FILE`` with a
two-column wavelength/density table. Only one telescope is supported: all telescopes of a run
receive the same photons.

Direction convention (docs/sources/presets.md#directions-zenith-and-azimuth): the zenith and the
astronomical azimuth (clockwise from north, east is 90 deg) give the direction the photons come
from, the telescope pointing. The default pointing (zenith 0) looks straight up.
"""

from __future__ import annotations

import argparse

from obdeect.sources.presets import _common
from obdeect.sources.presets._beam import non_negative, positive, unit_vector

NAME = "night-sky"
HELP = "night-sky background on one telescope sphere (sky_background emission per event)"


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the night-sky preset (telescope position, radius and view cone required)."""
    _common.add_run_arguments(parser)
    _common.add_profile_arguments(parser, angular=None)
    group = parser.add_argument_group("background")
    group.add_argument(
        "--photons", type=float, required=True, metavar="N", help="photons per event"
    )
    group.add_argument(
        "--bunch-size", type=float, default=1.0, metavar="F", help="photons per bunch (>= 1)"
    )
    group.add_argument(
        "--viewcone-deg",
        type=float,
        required=True,
        metavar="DEG",
        help="half opening angle of the viewed cone, in (0, 180] [deg]",
    )
    group = parser.add_argument_group("telescope")
    group.add_argument(
        "--telescope-position",
        type=_common.vector,
        required=True,
        metavar="X,Y,Z",
        help="telescope position in the sky frame [m]",
    )
    group.add_argument(
        "--telescope-radius",
        type=float,
        required=True,
        metavar="M",
        help="telescope sphere radius [m]",
    )
    group.add_argument(
        "--telescope-zenith",
        type=float,
        default=0.0,
        metavar="DEG",
        help="zenith of the sky direction the photons come from, the telescope pointing "
        "(default 0: up) [deg]",
    )
    group.add_argument(
        "--telescope-azimuth",
        type=float,
        default=0.0,
        metavar="DEG",
        help="astronomical azimuth of the sky direction the photons come from, the telescope "
        "pointing: clockwise from north, east is 90 [deg] (default 0)",
    )


def build(args: argparse.Namespace) -> dict:
    """Run description with one sky_background emission per event."""
    radius = positive(args.telescope_radius, "telescope radius")
    viewcone = float(args.viewcone_deg)
    if not (0.0 < viewcone <= 180.0):
        raise ValueError(f"view cone {viewcone} deg must lie in (0, 180]")
    photons = non_negative(args.photons, "photons")
    bunch = positive(args.bunch_size, "bunch size")
    if bunch < 1.0:
        raise ValueError(f"bunch size must be >= 1 (got {args.bunch_size})")

    pointing = unit_vector(args.telescope_zenith, args.telescope_azimuth)
    centre = list(args.telescope_position)
    emission = {
        "kind": "sky_background",
        "centre_m": centre,
        "radius_m": radius,
        "axis": [0.0 - x for x in pointing],
        "viewcone_deg": viewcone,
        "photons": photons,
        "bunch_size": bunch,
        "spectrum": _common.profile(args.spectrum, "spectrum"),
        "pulse": _common.profile(args.pulse, "pulse"),
    }
    return {
        "run": _common.run_section(args, frame="sky"),
        "telescopes": _common.single_telescope(radius, centre),
        "events": [{"repeat": args.events, "emissions": [emission]}],
    }
