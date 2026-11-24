"""Flasher: a light source in front of the camera (dish) or of a dual-mirror telescope.

Dish: one flasher at ``--position`` in the camera frame, pointing along -z, emitting into the
camera plane (``to_plane_z0``). Secondary: four flashers in the telescope frame; flasher ``i``
sits at ``(+-offset, +-offset, mirror_camera_distance + depth)`` and points at altitude
``-90 deg + inclination``, towards its azimuth.
"""

import argparse
import math

from obdeect.sources.presets import _common

NAME = "flasher"
HELP = "camera flasher: dish in front of the camera, or four flashers of a dual-mirror telescope"

_DOWN = [0.0, 0.0, -1.0]
_SECONDARY_ONLY = (
    "flasher_offset",
    "flasher_depth",
    "flasher_inclination",
    "mirror_camera_distance",
    "fiducial_radius",
    "delays",
    "fire",
)
_DISH_ONLY = ("position", "camera_radius")


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the flasher preset."""
    parser.add_argument("--type", choices=["dish", "secondary"], default="dish")
    _common.add_run_arguments(parser)
    _common.add_profile_arguments(parser, angular="isotropic")
    parser.add_argument(
        "--photons",
        type=float,
        nargs="+",
        required=True,
        metavar="N",
        help="photons into 4 pi per flasher; one event per value",
    )
    parser.add_argument(
        "--bunch-size", type=float, default=1.0, metavar="N", help="photons per bunch (>= 1)"
    )
    dish = parser.add_argument_group("dish (--type dish)")
    dish.add_argument(
        "--position",
        type=_common.vector,
        metavar="X,Y,Z",
        help="flasher position in the camera frame [m]",
    )
    dish.add_argument(
        "--camera-radius",
        type=float,
        metavar="M",
        help="camera radius [m]; sets the telescope radius of the target rule",
    )
    sec = parser.add_argument_group("secondary (--type secondary)")
    sec.add_argument(
        "--flasher-offset",
        type=float,
        metavar="M",
        help="flasher distance from the axis, in x and y [m]",
    )
    sec.add_argument(
        "--flasher-depth", type=float, metavar="M", help="flasher depth below the mirror plane [m]"
    )
    sec.add_argument(
        "--flasher-inclination",
        type=float,
        metavar="DEG",
        help="pointing is altitude -90 deg + inclination",
    )
    sec.add_argument(
        "--mirror-camera-distance",
        type=float,
        metavar="M",
        help="distance from the secondary mirror to the camera [m]",
    )
    sec.add_argument(
        "--fiducial-radius",
        type=float,
        metavar="M",
        help="telescope fiducial radius [m]; sets the target rule",
    )
    sec.add_argument(
        "--delays",
        type=float,
        nargs=4,
        metavar="NS",
        help="time offsets of flashers 1 to 4 [ns] (default 0)",
    )
    sec.add_argument(
        "--fire",
        choices=["all", "alternate", "1", "2", "3", "4"],
        help="all: one event with all flashers; alternate: event k fires "
        "flasher k mod 4; 1..4: only that flasher (default all)",
    )


def build(args: argparse.Namespace) -> dict:
    """Run description of the flasher."""
    if args.type == "dish":
        return _dish(args)
    return _secondary(args)


def _required(value, option: str):
    if value is None:
        raise ValueError(f"{option} is required")
    return value


def _reject(args: argparse.Namespace, names: tuple, applies: str) -> None:
    given = [f"--{name.replace('_', '-')}" for name in names if getattr(args, name) is not None]
    if given:
        raise ValueError(f"{', '.join(given)} applies only to --type {applies}")


def _profiles(args: argparse.Namespace) -> dict:
    return {
        "spectrum": _common.profile(args.spectrum, "spectrum"),
        "pulse": _common.profile(args.pulse, "pulse"),
        "angular": _common.profile(args.angular, "angular"),
    }


def _point(args, profiles, position, axis, photons, *, offset_ns=0.0, to_plane_z0=False):
    return _common.point_emission(
        photons=photons,
        bunch_size=args.bunch_size,
        position_m=position,
        axis=axis,
        angular=profiles["angular"],
        spectrum=profiles["spectrum"],
        pulse=profiles["pulse"],
        target="auto",
        time_offset_ns=offset_ns,
        to_plane_z0=to_plane_z0,
    )


def _dish(args: argparse.Namespace) -> dict:
    _reject(args, _SECONDARY_ONLY, "secondary")
    position = _required(args.position, "--position X,Y,Z (camera frame) with --type dish")
    radius = _required(args.camera_radius, "--camera-radius with --type dish")
    profiles = _profiles(args)
    events = [
        {
            "repeat": args.events,
            "emissions": [_point(args, profiles, position, _DOWN, photons, to_plane_z0=True)],
        }
        for photons in args.photons
    ]
    return {
        "run": _common.run_section(args, frame="camera"),
        "telescopes": _common.single_telescope(radius),
        "events": events,
    }


def _flasher(index: int, offset: float, depth: float, inclination: float, mirror: float):
    """Position and pointing of flasher ``index`` (0..3) in the telescope frame."""
    sx = 1.0 if index < 2 else -1.0
    sy = 1.0 if index % 2 == 0 else -1.0
    x, y = sx * offset, sy * offset
    azimuth = math.atan2(y, x)
    altitude = math.radians(-90.0 + inclination)
    axis = [
        math.cos(altitude) * math.cos(azimuth),
        math.cos(altitude) * math.sin(azimuth),
        math.sin(altitude),
    ]
    return [x, y, mirror + depth], axis


def _secondary(args: argparse.Namespace) -> dict:
    _reject(args, _DISH_ONLY, "dish")
    offset = _required(args.flasher_offset, "--flasher-offset with --type secondary")
    depth = _required(args.flasher_depth, "--flasher-depth with --type secondary")
    inclination = _required(args.flasher_inclination, "--flasher-inclination with --type secondary")
    mirror = _required(
        args.mirror_camera_distance, "--mirror-camera-distance with --type secondary"
    )
    radius = _required(args.fiducial_radius, "--fiducial-radius with --type secondary")
    delays = args.delays or [0.0, 0.0, 0.0, 0.0]
    fire = args.fire or "all"
    profiles = _profiles(args)
    geometry = [_flasher(i, offset, depth, inclination, mirror) for i in range(4)]

    def emission(index: int, photons: float) -> dict:
        position, axis = geometry[index]
        return _point(args, profiles, position, axis, photons, offset_ns=delays[index])

    events = []
    for photons in args.photons:
        if fire == "alternate":
            for k in range(args.events):
                index = k % 4
                events.append({
                    "repeat": 1,
                    "metadata": {"flashers": [float(index + 1)]},
                    "emissions": [emission(index, photons)],
                })
        else:
            fired = list(range(4)) if fire == "all" else [int(fire) - 1]
            events.append({
                "repeat": args.events,
                "metadata": {"flashers": [float(i + 1) for i in fired]},
                "emissions": [emission(i, photons) for i in fired],
            })
    return {
        "run": _common.run_section(args, frame="telescope"),
        "telescopes": _common.single_telescope(radius),
        "events": events,
    }
