"""Options and helpers shared by the presets.

Profile options take ``KIND`` or ``KIND:NAME=VALUE,...`` with the keys of the run description,
and ``table:FILE`` for a two-column table file::

    --spectrum fixed:wavelength_nm=400
    --spectrum inverse_square:low_nm=300,high_nm=600
    --pulse gauss:sigma_ns=2.5
    --angular cone:angle_deg=1.2
    --spectrum table:night_sky.dat
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

#: Philox stage of the random choices of presets (obdeect stream ``stage_programs``).
STAGE_PRESETS = 0x444C0002

# Keys of the two table columns per profile.
_TABLE_KEYS = {
    "spectrum": ("wavelength_nm", "density"),
    "pulse": ("time_ns", "rate"),
    "angular": ("angle_deg", "value"),
}


def read_table(path: str | Path) -> tuple[list[float], list[float]]:
    """Two-column table: ``#`` starts a comment, blank lines are skipped, others need 2 numbers."""
    x, y = [], []
    for number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        text = line.split("#", 1)[0].split()
        if not text:
            continue
        if len(text) < 2:
            raise ValueError(f"{path} line {number}: needs two numbers")
        try:
            x.append(float(text[0]))
            y.append(float(text[1]))
        except ValueError:
            raise ValueError(f"{path} line {number}: not a number") from None
    if not x:
        raise ValueError(f"{path}: no rows")
    return x, y


def profile(text: str, section: str) -> dict:
    """Profile dictionary of a ``KIND[:NAME=VALUE,...]`` or ``table:FILE`` option."""
    kind, _, rest = text.partition(":")
    if not kind:
        raise ValueError(f"{section} '{text}': needs a kind")
    if kind == "table":
        if not rest:
            raise ValueError(f"{section} '{text}': table needs a file")
        x_key, y_key = _TABLE_KEYS[section]
        x, y = read_table(rest)
        return {"kind": "table", x_key: x, y_key: y}
    out: dict = {"kind": kind}
    for item in filter(None, rest.split(",")):
        name, eq, value = item.partition("=")
        if not eq or not name:
            raise ValueError(f"{section} '{text}': '{item}' is not NAME=VALUE")
        try:
            out[name] = float(value)
        except ValueError:
            raise ValueError(f"{section} '{text}': '{value}' is not a number") from None
    return out


def vector(text: str) -> list[float]:
    """``x,y,z`` as three floats."""
    items = text.split(",")
    if len(items) != 3:
        raise argparse.ArgumentTypeError(f"'{text}' needs three comma-separated values")
    try:
        return [float(v) for v in items]
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{text}' is not three numbers") from None


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the run section and the event count."""
    group = parser.add_argument_group("run")
    group.add_argument(
        "--observation-level",
        type=float,
        required=True,
        metavar="M",
        help="observation level above sea level [m]",
    )
    group.add_argument("--run-number", type=int, default=1, metavar="N")
    group.add_argument("--seed", type=int, default=1, metavar="N")
    group.add_argument(
        "--atmosphere-id",
        type=int,
        default=-1,
        metavar="N",
        help="atmosphere profile number (-1: none)",
    )
    group.add_argument("--events", type=int, default=1, metavar="N", help="events to write")


def add_profile_arguments(parser: argparse.ArgumentParser, *, angular: str | None) -> None:
    """--spectrum (required), --pulse, and --angular when the preset takes one."""
    group = parser.add_argument_group("profiles")
    group.add_argument("--spectrum", required=True, metavar="KIND[:K=V,...]")
    group.add_argument("--pulse", default="fixed:time_ns=0", metavar="KIND[:K=V,...]")
    if angular is not None:
        group.add_argument("--angular", default=angular, metavar="KIND[:K=V,...]")


def run_section(args: argparse.Namespace, frame: str = "sky", **metadata) -> dict:
    """The run section from the run options."""
    run = {
        "number": args.run_number,
        "seed": args.seed,
        "observation_level_m": args.observation_level,
        "atmosphere_id": args.atmosphere_id,
        "frame": frame,
    }
    if metadata:
        run["metadata"] = metadata
    return run


def choice_rng(seed: int) -> np.random.Generator:
    """Generator of a preset's random choices (pixels, core positions), never of photons.

    One rule for all presets: Philox keyed on (run seed, ``STAGE_PRESETS``); photons use the
    emission seeds of the run (``source_seed(run.seed, k)``), so the two never share draws.
    """
    return np.random.Generator(np.random.Philox(key=[seed % 2**64, STAGE_PRESETS]))
