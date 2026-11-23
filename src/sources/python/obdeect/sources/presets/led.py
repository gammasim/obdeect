"""LED: light sources in front of selected camera pixels (camera frame).

Pixels come from a camera definition file: ``Pixel <id> <type> <x> <y>`` lines with positions in
cm (the file's own unit), and ``Trigger`` lines naming groups of pixels. Pixel lists look like
``1 +2 3 [4,5]``: ``+`` marks a required pixel, ``[..]`` lists pre-summed slave pixels, ``%``
starts a comment. Each lit pixel is one point emission, parallel along -z, into the camera plane.
Old program: ``led.hpp``.
"""

import argparse
import re
import warnings
from itertools import combinations
from pathlib import Path
from typing import NamedTuple

from obdeect.sources.presets import _common

NAME = "led"
HELP = "LEDs in front of selected camera pixels (camera frame)"

_PARALLEL = {"kind": "cone", "angle_deg": 0.0}
_ITEM = re.compile(r"([+]?)[ \t]*(\d+)[ \t]*(?:\[([^\]]*)\])?")


class _Item(NamedTuple):
    id: int
    required: bool
    slaves: list[int]


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the LED preset."""
    parser.add_argument(
        "--camera-file",
        required=True,
        metavar="FILE",
        help="camera definition: Pixel and Trigger lines (positions in cm)",
    )
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--use", metavar="LIST", help="pixel list, e.g. '1 +2 [3,4] 5'")
    modes.add_argument(
        "--use-trg-rand", type=int, metavar="N", help="N events from random trigger groups"
    )
    modes.add_argument(
        "--use-trg-all", action="store_true", help="one list per trigger group of the camera file"
    )
    parser.add_argument(
        "--require",
        type=int,
        default=0,
        metavar="N",
        help="pixels per event of each list (0: all pixels of the list)",
    )
    parser.add_argument(
        "--random", action="store_true", help="pick the --require pixels of each event at random"
    )
    parser.add_argument(
        "--slaves", action="store_true", help="also fire the pre-summed slave pixels"
    )
    _common.add_run_arguments(parser)
    _common.add_profile_arguments(parser, angular=None)
    parser.add_argument(
        "--photons", type=float, required=True, metavar="N", help="photons into 4 pi per LED"
    )
    parser.add_argument(
        "--bunch-size", type=float, default=1.0, metavar="N", help="photons per bunch (>= 1)"
    )
    parser.add_argument(
        "--camera-radius",
        type=float,
        required=True,
        metavar="M",
        help="camera radius [m]; sets the telescope radius of the target rule",
    )


def build(args: argparse.Namespace) -> dict:
    """Run description: one event per pixel combination."""
    if args.random and args.require <= 0:
        raise ValueError("--random needs --require N with N > 0")
    pixels, triggers = read_camera(args.camera_file, args.slaves)
    rng = _common.choice_rng(args.seed)
    spectrum = _common.profile(args.spectrum, "spectrum")
    pulse = _common.profile(args.pulse, "pulse")
    events = []
    for items in _pixel_lists(args, triggers, rng):
        _check_ids(items, pixels, args.slaves)
        size = len(items)
        nreq = args.require or size
        if nreq > size:
            warnings.warn(f"required {nreq} pixels but the list has {size}; list skipped")
            continue
        if nreq == size:
            events.append(_event(args, items, pixels, spectrum, pulse, args.events))
        elif args.random:
            for _ in range(args.events):
                picks = rng.choice(size, size=nreq, replace=False)
                events.append(_event(args, [items[k] for k in picks], pixels, spectrum, pulse, 1))
        else:
            for combo in _combinations(items, nreq):
                events.append(_event(args, combo, pixels, spectrum, pulse, args.events))
    return {
        "run": _common.run_section(args, frame="camera"),
        "telescopes": [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": args.camera_radius}],
        "events": events,
    }


def read_camera(path: str, keep_slaves: bool) -> tuple[dict[int, tuple[float, float]], list[str]]:
    """Pixel positions [cm] by id, and the pixel-list text of each trigger group."""
    pixels: dict[int, tuple[float, float]] = {}
    triggers: list[str] = []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        words = line.split()
        keyword = words[0].lower()
        if keyword == "pixel":
            try:
                pixel_id, x, y = int(words[1]), float(words[3]), float(words[4])
            except IndexError, ValueError:
                warnings.warn(f"problem with pixel definition, line skipped: {line}")
                continue
            if pixel_id < 0:
                raise ValueError(f"negative pixel id {pixel_id}")
            if pixel_id in pixels:
                raise ValueError(f"duplicate pixel id {pixel_id}")
            pixels[pixel_id] = (x, y)
        elif keyword in ("rotate", "pixtype"):
            continue
        elif "Trigger" in line:
            triggers.append(trigger_text(line, keep_slaves))
    return pixels, triggers


def trigger_text(line: str, keep_slaves: bool) -> str:
    """Pixel-list text of a trigger line: the text after the first 'of', without slave groups."""
    text = line
    position = text.find("of")
    if position >= 0:
        text = text[position + 2 :].lstrip(" \t")
    return text if keep_slaves else re.sub(r"\[[^\]]*\]", "", text)


def parse_pixel_list(text: str) -> list[_Item]:
    """Items of a pixel list; a repeated pixel is used once, with a warning."""
    rest = re.split(r"[%\n]", text, maxsplit=1)[0].strip(" \t")
    items: list[_Item] = []
    seen: set[int] = set()
    while rest:
        match = _ITEM.match(rest)
        if match is None:
            raise ValueError(f"invalid in list of pixels to use: '{rest}'")
        slaves = match.group(3)
        if slaves is not None and not re.fullmatch(r"[\d \t,]*", slaves):
            raise ValueError(f"invalid pre-sum pixel group [{slaves}] in '{text}'")
        rest = rest[match.end() :].lstrip(" \t")
        item = _Item(
            int(match.group(2)),
            bool(match.group(1)),
            [int(s) for s in re.findall(r"\d+", slaves or "")],
        )
        if item.id in seen:
            warnings.warn(f"duplicated pixel {item.id} in use list, used only once")
            continue
        seen.add(item.id)
        items.append(item)
    if not items:
        raise ValueError("failed to obtain a pixel list")
    return items


def _pixel_lists(args, triggers: list[str], rng) -> list[list[_Item]]:
    if args.use is not None:
        return [parse_pixel_list(trigger_text(args.use, args.slaves))]
    if not triggers:
        raise ValueError("the camera definition has no trigger groups")
    if args.use_trg_all:
        return [parse_pixel_list(text) for text in triggers]
    return [
        parse_pixel_list(triggers[int(rng.integers(len(triggers)))])
        for _ in range(args.use_trg_rand)
    ]


def _check_ids(items: list[_Item], pixels: dict, slaves_kept: bool) -> None:
    for item in items:
        ids = [item.id, *(item.slaves if slaves_kept else [])]
        for pixel_id in ids:
            if pixel_id not in pixels:
                raise ValueError(f"pixel ID {pixel_id} is not defined in the camera file")


def _combinations(items: list[_Item], nreq: int) -> list[list[_Item]]:
    """Combinations of nreq pixels: the required ones first, then optional ones in list order."""
    required = [k for k, item in enumerate(items) if item.required]
    optional = [k for k, item in enumerate(items) if not item.required]
    if len(required) >= nreq:
        picks = [required]
    else:
        picks = [
            sorted([*required, *extra]) for extra in combinations(optional, nreq - len(required))
        ]
    return [[items[k] for k in pick] for pick in picks]


def _event(args, items: list[_Item], pixels: dict, spectrum, pulse, repeat: int) -> dict:
    """One event firing the given pixels: a point emission per pixel (and per slave, if kept)."""
    emissions = []
    for item in items:
        for pixel_id in [item.id, *(item.slaves if args.slaves else [])]:
            x, y = pixels[pixel_id]
            emissions.append({
                "kind": "point",
                "photons": args.photons,
                "bunch_size": args.bunch_size,
                "time_offset_ns": 0.0,
                "to_plane_z0": False,
                "position_m": [x / 100.0, y / 100.0, 0.0],
                "axis": [0.0, 0.0, -1.0],
                "angular": _PARALLEL,
                "spectrum": spectrum,
                "pulse": pulse,
                "target": "auto",
            })
    labels = [float(item.id) for item in items]
    return {"repeat": repeat, "metadata": {"pixels": labels}, "emissions": emissions}
