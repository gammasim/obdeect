"""Convert older extinction tables to the obdeect propagation transmission-table format.

Both legacy layouts give a cumulative vertical optical depth c(lambda, h) that rises with
altitude. The optical depth between an observation level h_obs and an emission height h_e is
c(h_e) - c(h_obs), with c linearly interpolated in height and in wavelength. The converter
writes tau(lambda, h) = c(lambda, h) - c(lambda, h_obs) for integer wavelengths (1 nm steps)
on a height grid above h_obs; the absolute offset of c does not matter. See
docs/extinction.md, "Legacy tables".
"""

import argparse
import math
import os
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

STEP_NM = 5.0
MAX_WAVELENGTH_NM = 999

#: Names accepted for ``format``, with the default minimum wavelength [nm] (None: in the file).
FORMATS: dict[str, float | None] = {
    "corsika-atmabs": None,
    "kascade": 180.0,
    "modtran4": 200.0,
    "modtran5": 205.0,
}

MAX_HEIGHTS = 100  # most heights of a obdeect propagation table, H2 included
NEAR_OFFSETS_KM = (0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)

# corsika-atmabs: heights 0, 1, ..., 50 km per wavelength.
ATMABS_HEIGHTS_KM = 51


@dataclass(frozen=True)
class LegacyTable:
    """Cumulative vertical optical depth ``depth[i_wavelength, i_altitude]``."""

    wavelength_nm: np.ndarray
    altitude_km: np.ndarray
    depth: np.ndarray


def _number(token: str, where: str) -> float:
    try:
        value = float(token)
    except ValueError:
        raise ValueError(f"{where}: {token!r} is not a number") from None
    if not math.isfinite(value):
        raise ValueError(f"{where}: {token!r} is not finite")
    return value


def _data_lines(path: str | os.PathLike, skip_extinction: bool):
    with open(path, encoding="utf-8") as stream:
        for number, line in enumerate(stream, start=1):
            text = line.split("#", 1)[0].strip()
            if not text or (skip_extinction and "EXTINCTION" in text):
                continue
            yield number, text.split()


def _check_uniform(wavelength_nm: np.ndarray, name: str) -> None:
    if len(wavelength_nm) < 2 or not np.allclose(
        np.diff(wavelength_nm), STEP_NM, rtol=0.0, atol=1e-9
    ):
        raise ValueError(f"{name}: wavelengths must rise in steps of {STEP_NM:g} nm, at least two")


def read_corsika_atmabs(path: str | os.PathLike) -> LegacyTable:
    """Blocks of one wavelength [nm] line and 51 values (0 to 50 km), 'EXTINCTION' header lines."""
    wavelengths: list[float] = []
    rows: list[list[float]] = []
    pending: list[float] | None = None
    for number, tokens in _data_lines(path, skip_extinction=True):
        where = f"{path} line {number}"
        if pending is None or len(pending) == ATMABS_HEIGHTS_KM:
            if pending is not None:
                rows.append(pending)
            if len(tokens) != 1:
                raise ValueError(f"{where}: expected a line with one wavelength")
            wavelengths.append(_number(tokens[0], where))
            pending = []
            continue
        pending.extend(_number(token, where) for token in tokens)
        if len(pending) > ATMABS_HEIGHTS_KM:
            raise ValueError(f"{where}: more than {ATMABS_HEIGHTS_KM} values for one wavelength")
    if pending is None:
        raise ValueError(f"{path}: no data")
    if len(pending) != ATMABS_HEIGHTS_KM:
        raise ValueError(f"{path}: last block has {len(pending)} of {ATMABS_HEIGHTS_KM} values")
    rows.append(pending)
    wavelength_nm = np.array(wavelengths)
    _check_uniform(wavelength_nm, str(path))
    return LegacyTable(wavelength_nm, np.arange(ATMABS_HEIGHTS_KM, dtype=float), np.array(rows))


def read_matrix(path: str | os.PathLike, min_wavelength_nm: float) -> LegacyTable:
    """One row per 1 km step from 0 km, one column per 5 nm step from ``min_wavelength_nm``."""
    rows: list[list[float]] = []
    for number, tokens in _data_lines(path, skip_extinction=False):
        where = f"{path} line {number}"
        rows.append([_number(token, where) for token in tokens])
        if len(rows[-1]) != len(rows[0]):
            raise ValueError(f"{where}: {len(rows[-1])} columns, expected {len(rows[0])}")
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError(f"{path}: needs at least two altitude rows and two wavelength columns")
    depth = np.array(rows).T  # [wavelength, altitude]
    wavelength_nm = min_wavelength_nm + STEP_NM * np.arange(depth.shape[0])
    return LegacyTable(wavelength_nm, np.arange(depth.shape[1], dtype=float), depth)


def read_legacy_table(
    path: str | os.PathLike, format: str, min_wavelength_nm: float | None = None
) -> LegacyTable:
    """Read a legacy table; ``format`` is one of :data:`FORMATS`."""
    if format not in FORMATS:
        raise ValueError(f"unknown format {format!r}; use one of {', '.join(FORMATS)}")
    if format == "corsika-atmabs":
        if min_wavelength_nm is not None:
            raise ValueError("corsika-atmabs files hold their wavelengths; no minimum wavelength")
        return read_corsika_atmabs(path)
    minimum = FORMATS[format] if min_wavelength_nm is None else min_wavelength_nm
    if not (isinstance(minimum, float | int) and math.isfinite(minimum) and minimum > 0):
        raise ValueError("minimum wavelength must be a positive finite number")
    return read_matrix(path, float(minimum))


def height_grid_km(altitude_km: np.ndarray, observation_km: float) -> np.ndarray:
    """Heights above the observation level: every source altitude, eighth steps in the lowest
    6 km above it, and a few points within 500 m of it (tau is linear in log10 of the height).
    Coarser steps (4, 2, 1 parts) are used if the grid would exceed ``MAX_HEIGHTS``."""
    knots = [observation_km, *(a for a in altitude_km if a > observation_km)]
    for parts_low in (8, 4, 2, 1):
        grid: list[float] = []
        for low, high in zip(knots[:-1], knots[1:], strict=True):
            parts = parts_low if low < observation_km + 6.0 else 1
            grid.extend(low + (high - low) * k / parts for k in range(1, parts + 1))
        near = [observation_km + offset for offset in NEAR_OFFSETS_KM if parts_low > 1]
        grid = [*(h for h in near if h < grid[0]), *grid]
        if len(grid) < MAX_HEIGHTS:  # one slot is H2
            return np.array(grid)
    raise ValueError(
        f"the table has {len(knots) - 1} altitudes above the observation level; "
        f"at most {MAX_HEIGHTS - 1} heights above it are supported"
    )


def convert_table(
    input_path: str | os.PathLike,
    output_path: str | os.PathLike,
    *,
    format: str,
    observation_level_m: float,
    min_wavelength_nm: float | None = None,
) -> Path:
    """Convert a legacy extinction table to the obdeect propagation transmission-table format.

    Args:
        input_path: legacy table.
        output_path: obdeect propagation table (``# H2= ..., H1= ...``) to write.
        format: ``corsika-atmabs``, ``kascade``, ``modtran4`` or ``modtran5``.
        observation_level_m: observation level [m a.s.l.] H2; the optical depth is counted from
            it. Must lie inside the altitude range of the table.
        min_wavelength_nm: first wavelength of the matrix formats (default: 180, 200, 205 nm for
            kascade, modtran4, modtran5).

    Returns:
        The output path. Wavelengths are integer nm from the first to the last source
        wavelength (at most 999 nm), optical depths are linear in wavelength between the source
        columns.
    """
    if not math.isfinite(observation_level_m):
        raise ValueError("observation level must be finite")
    table = read_legacy_table(input_path, format, min_wavelength_nm)
    observation_km = observation_level_m / 1000.0
    if not (table.altitude_km[0] <= observation_km < table.altitude_km[-1]):
        raise ValueError(
            f"observation level {observation_level_m:g} m is outside the table altitudes "
            f"{table.altitude_km[0]:g} to {table.altitude_km[-1]:g} km"
        )
    if observation_km <= 0.0:
        raise ValueError("observation level must be above sea level (H2 > 0)")
    heights = height_grid_km(table.altitude_km, observation_km)
    wavelengths = np.arange(
        math.ceil(table.wavelength_nm[0]),
        min(math.floor(table.wavelength_nm[-1]), MAX_WAVELENGTH_NM) + 1,
    )
    if len(wavelengths) == 0:
        raise ValueError("no integer wavelength below 1000 nm in the table")
    # Linear in height for every source wavelength, then linear in wavelength.
    at_heights = np.array([np.interp(heights, table.altitude_km, row) for row in table.depth])
    at_observation = np.array([
        np.interp(observation_km, table.altitude_km, row) for row in table.depth
    ])
    tau_source = at_heights - at_observation[:, None]
    above = table.altitude_km > observation_km
    chain = np.column_stack([at_observation, table.depth[:, above]])
    tolerance = 1e-9 * np.maximum(1.0, np.maximum(abs(chain[:, :-1]), abs(chain[:, 1:])))
    if (chain[:, 1:] < chain[:, :-1] - tolerance).any():
        raise ValueError("optical depth decreases with altitude above the observation level")
    tau_source = np.clip(tau_source, 0.0, None)
    lines = [
        "# Converted by obdeect.propagation.convert_table",
        f"# Source format: {format}",
        "# tau(lambda, h) = c(lambda, h) - c(lambda, H2), vertical optical depth above H2",
        f"# H2= {observation_km:.6f}, H1= " + " ".join(f"{h:.6f}" for h in heights),
    ]
    for wavelength in wavelengths:
        values = [
            float(np.interp(wavelength, table.wavelength_nm, tau_source[:, k]))
            for k in range(len(heights))
        ]
        lines.append(f"{wavelength:d} " + " ".join(f"{v:.12g}" for v in values))
    output = Path(output_path)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def main(argv: Sequence[str] | None = None) -> int:
    """Command line: ``obdeect propagation convert-table --format FORMAT IN -o OUT
    --observation-level M``."""
    parser = argparse.ArgumentParser(
        prog="obdeect propagation convert-table",
        description=(
            "Convert a legacy extinction table to the obdeect propagation "
            "transmission-table format."
        ),
    )
    parser.add_argument("input", help="legacy table")
    parser.add_argument("-o", "--output", required=True, help="obdeect propagation table to write")
    parser.add_argument("--format", required=True, choices=list(FORMATS))
    parser.add_argument(
        "--observation-level",
        type=float,
        required=True,
        metavar="M",
        help="observation level [m a.s.l.]",
    )
    parser.add_argument(
        "--min-wavelength",
        type=float,
        metavar="NM",
        help="first wavelength of the matrix formats (default 180 kascade, 200 modtran4, "
        "205 modtran5)",
    )
    args = parser.parse_args(argv)
    try:
        convert_table(
            args.input,
            args.output,
            format=args.format,
            observation_level_m=args.observation_level,
            min_wavelength_nm=args.min_wavelength,
        )
    except (OSError, ValueError) as error:
        print(f"obdeect propagation convert-table: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
