"""Point source: a light source at a finite distance above the observation level (sky frame).

``--position`` is in the sky frame with z the height above the observation level. The emitted
direction is ``--direction`` normalised (default: downward). Each ``--photons`` value is one event.
Old program: ``point-source.hpp``.
"""

import argparse
import math

from obdeect.sources.presets import _common

NAME = "point-source"
HELP = "point light source at a finite distance above the observation level (sky frame)"


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the point-source preset."""
    _common.add_run_arguments(parser)
    _common.add_profile_arguments(parser, angular="isotropic")
    parser.add_argument(
        "--position",
        type=_common.vector,
        required=True,
        metavar="X,Y,Z",
        help="source position in the sky frame [m]; z is the height above the observation level",
    )
    parser.add_argument(
        "--direction",
        type=_common.vector,
        default=[0.0, 0.0, -1.0],
        metavar="X,Y,Z",
        help="emission direction, normalised (default: down)",
    )
    parser.add_argument(
        "--photons",
        type=float,
        nargs="+",
        required=True,
        metavar="N",
        help="photons into 4 pi per event; one event per value",
    )
    parser.add_argument(
        "--bunch-size", type=float, default=1.0, metavar="N", help="photons per bunch (>= 1)"
    )
    where = parser.add_mutually_exclusive_group(required=True)
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


def _telescope(text: str) -> dict:
    items = text.split(",")
    try:
        x, y, z, r = (float(v) for v in items)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{text}' needs four numbers X,Y,Z,R") from None
    return {"position_m": [x, y, z], "radius_m": r}


def build(args: argparse.Namespace) -> dict:
    """Run description: one point emission per ``--photons`` value."""
    norm = math.sqrt(sum(v * v for v in args.direction))
    if norm == 0.0:
        raise ValueError("--direction must not be the zero vector")
    axis = [v / norm for v in args.direction]
    spectrum = _common.profile(args.spectrum, "spectrum")
    pulse = _common.profile(args.pulse, "pulse")
    angular = _common.profile(args.angular, "angular")
    events = [
        {
            "repeat": args.events,
            "emissions": [
                {
                    "kind": "point",
                    "photons": photons,
                    "bunch_size": args.bunch_size,
                    "time_offset_ns": 0.0,
                    "to_plane_z0": False,
                    "position_m": list(args.position),
                    "axis": axis,
                    "angular": angular,
                    "spectrum": spectrum,
                    "pulse": pulse,
                    "target": "auto",
                }
            ],
        }
        for photons in args.photons
    ]
    return {
        "run": _common.run_section(args, frame="sky"),
        "telescopes": [
            {"id": i, **t}
            for i, t in enumerate(
                args.telescope
                or [{"position_m": [0.0, 0.0, 0.0], "radius_m": args.telescope_radius}],
                start=1,
            )
        ],
        "events": events,
    }
