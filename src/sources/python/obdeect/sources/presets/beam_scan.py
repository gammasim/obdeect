"""Beam scan: a parallel beam moved along two lines over the camera area (sky frame).

The scan starts at ``--scan-start``, runs along x for ``--scan-extent-x`` and then along y for
``--scan-extent-y``, in steps of ``--scan-step``. Each scan position is one point emission of a
single event, at ``--height``. The beam points downwards by default; ``--beam-zenith`` (zenith,
> 90 down, < 90 up) and ``--beam-azimuth`` set another direction. A horizontal beam is rejected:
it never reaches the camera plane (planned with the sphere representation). Old program:
``beam-scan.hpp`` (its extents and step were fractions of the radius).
"""

import argparse
import math

from obdeect.sources.presets import _common
from obdeect.sources.presets._beam import unit_vector

NAME = "beam-scan"
HELP = "parallel beam scanned along two lines over the camera area (sky frame)"

_PARALLEL = {"kind": "cone", "angle_deg": 0.0}
_DOWNWARD_ZENITH = 180.0


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the beam-scan preset."""
    _common.add_run_arguments(parser)
    _common.add_profile_arguments(parser, angular=None)
    parser.add_argument(
        "--photons-per-position",
        type=float,
        required=True,
        metavar="N",
        help="photons into 4 pi per scan position",
    )
    parser.add_argument(
        "--bunch-size", type=float, default=1.0, metavar="N", help="photons per bunch (>= 1)"
    )
    parser.add_argument(
        "--camera-radius",
        type=float,
        required=True,
        metavar="M",
        help="camera radius [m]; telescope radius at the origin",
    )
    beam = parser.add_argument_group("beam")
    beam.add_argument(
        "--beam-zenith",
        type=float,
        default=_DOWNWARD_ZENITH,
        metavar="DEG",
        help="zenith of the beam: 180 (default) points down, 0 up; 90 (horizontal) is rejected "
        "[deg]",
    )
    beam.add_argument(
        "--beam-azimuth",
        type=float,
        default=0.0,
        metavar="DEG",
        help="astronomical azimuth of the beam: clockwise from north (default 0) [deg]",
    )
    scan = parser.add_argument_group("scan")
    scan.add_argument(
        "--scan-start",
        type=float,
        nargs=2,
        required=True,
        metavar=("X", "Y"),
        help="first scan position in the sky frame [m]",
    )
    scan.add_argument(
        "--scan-extent-x",
        type=float,
        required=True,
        metavar="M",
        help="length of the scan line along x [m]",
    )
    scan.add_argument(
        "--scan-extent-y",
        type=float,
        required=True,
        metavar="M",
        help="length of the scan line along y [m]",
    )
    scan.add_argument(
        "--scan-step",
        type=float,
        required=True,
        metavar="M",
        help="distance between scan positions [m]",
    )
    scan.add_argument(
        "--height",
        type=float,
        required=True,
        metavar="M",
        help="beam source height in the sky frame [m]",
    )


def scan_positions(
    start: tuple[float, float], extent_x: float, extent_y: float, step: float
) -> list[tuple[float, float]]:
    """Positions along x at the start y, then along y at the start x (start once)."""
    if not step > 0.0:
        raise ValueError("--scan-step must be > 0")
    if extent_x < 0.0 or extent_y < 0.0:
        raise ValueError("--scan-extent-x and --scan-extent-y must be >= 0")
    x0, y0 = start
    positions = [(x0 + k * step, y0) for k in range(math.ceil(extent_x / step))]
    first_y = 1 if positions else 0  # the start is already on the x line
    positions += [(x0, y0 + k * step) for k in range(first_y, math.ceil(extent_y / step))]
    if not positions:
        raise ValueError("the scan has no positions (check the extents and the start)")
    return positions


def build(args: argparse.Namespace) -> dict:
    """Run description: one event with a point emission per scan position."""
    axis = list(unit_vector(args.beam_zenith, args.beam_azimuth))
    if abs(axis[2]) <= 1e-12:
        raise ValueError(
            f"--beam-zenith {args.beam_zenith} deg is horizontal: horizontal beams are not "
            "supported yet (the photons never reach the camera plane)"
        )
    positions = scan_positions(
        tuple(args.scan_start), args.scan_extent_x, args.scan_extent_y, args.scan_step
    )
    spectrum = _common.profile(args.spectrum, "spectrum")
    pulse = _common.profile(args.pulse, "pulse")
    emissions = [
        {
            "kind": "point",
            "photons": args.photons_per_position,
            "bunch_size": args.bunch_size,
            "time_offset_ns": 0.0,
            "to_plane_z0": False,
            "position_m": [x, y, args.height],
            "axis": axis,
            "angular": _PARALLEL,
            "spectrum": spectrum,
            "pulse": pulse,
            "target": "none",
        }
        for x, y in positions
    ]
    return {
        "run": _common.run_section(args, frame="sky"),
        "telescopes": [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": args.camera_radius}],
        "events": [{"repeat": args.events, "emissions": emissions}],
    }
