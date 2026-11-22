"""lidar: a pulsed satellite lidar pointing down (or up) at the observation site.

Model: one point emission of the satellite's laser pulse, ``E lambda / (h c)`` primary photons
with pulse energy ``E`` (``--pulse-energy``) and wavelength ``--wavelength-nm``. The emission
starts at the satellite, or at ``--start-height`` on the beam (e.g. 100000 m: the atmosphere
above is negligible), with the beam axis, the divergence (full angle) and the pulse profile of
the options. Scattering is computed by propagation. Each pulse is one event (``--events``).

Satellite: ``--satellite-altitude M`` above the observation level, at the horizontal offset
``--satellite-offset X,Y`` (default 0,0), or ``--satellite-position X,Y,Z``.

Chain::

    obdeect sources run lidar ... -o primaries.pst
    obdeect propagation run --transmission T --site-altitude H --seed N --mode weighted \\
        --scattering S --scattering-points 64 -i primaries.pst -o scattered.pst

t = 0 is the laser firing. Limits: the beam is not horizontal and its tilt from the vertical
(nadir for a downward beam, the zenith for an upward one) is at most ``--max-tilt``; the
atmosphere is plane-parallel and only single scattering is computed.
Direction convention and the direct-beam check: see ``_pulsed_beam``.
"""

from __future__ import annotations

import argparse

from obdeect.sources.presets import _common, _pulsed_beam

NAME = "lidar"
HELP = "pulsed satellite lidar (up or down): primary photons of one pulse (point emission)"


def _offset(text: str) -> tuple[float, float]:
    items = text.split(",")
    try:
        x, y = (float(v) for v in items)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{text}' needs two numbers X,Y") from None
    return x, y


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the lidar preset (no device defaults: all geometry is required)."""
    _common.add_run_arguments(parser)
    group = parser.add_argument_group("satellite")
    where = group.add_mutually_exclusive_group(required=True)
    where.add_argument(
        "--satellite-position",
        type=_common.vector,
        metavar="X,Y,Z",
        help="satellite position in the sky frame; z above the observation level [m]",
    )
    where.add_argument(
        "--satellite-altitude",
        type=float,
        metavar="M",
        help="satellite height above the observation level [m]",
    )
    group.add_argument(
        "--satellite-offset",
        type=_offset,
        default=None,
        metavar="X,Y",
        help="horizontal offset of the satellite with --satellite-altitude (default 0,0) [m]",
    )
    _pulsed_beam.add_arguments(parser)


def build(args: argparse.Namespace) -> dict:
    """Run description of one pulse per event."""
    if args.satellite_position is not None:
        if args.satellite_offset is not None:
            raise ValueError("--satellite-offset is only used with --satellite-altitude")
        position = args.satellite_position
    else:
        x, y = args.satellite_offset or (0.0, 0.0)
        position = [x, y, args.satellite_altitude]
    return _pulsed_beam.build(args, position=position)
