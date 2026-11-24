"""sim_telarray conventions: RPOL tables, mirror lists, segmentations and structure tables.

Every reader converts to SI units and rejects input it cannot represent exactly; none of
them places optics in a telescope (see ``simulation_models``)."""

from __future__ import annotations

import bisect
import math
import re
from typing import Any

from obdeect.ray_tracing.errors import ModelError


class TableImportError(ValueError):
    """An RPOL table cannot be represented without discarding a declared option."""


def _polynomial_product(first, second):
    result = [0.0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return result


def _quadratic(x, values, origin, width):
    result = [0.0] * 3
    for i in range(3):
        polynomial = [values[i]]
        for j in range(3):
            if i != j:
                polynomial = _polynomial_product(
                    polynomial, [(origin - x[j]) / (x[i] - x[j]), width / (x[i] - x[j])]
                )
        for j, coefficient in enumerate(polynomial):
            result[j] += coefficient
    return result


def interval_coefficients(axis, values, scheme):
    """Normalized-interval cubic coefficients, resolved once before native tracing."""
    if scheme in (0, 1):
        return []
    count = len(axis)
    if count < 3:
        return [[a, b - a, 0.0, 0.0] for a, b in zip(values, values[1:])]
    if scheme == 2:
        coefficients = []
        for i in range(count - 1):
            left = max(0, i - 1)
            right = min(i, count - 3)
            width = axis[i + 1] - axis[i]
            a = _quadratic(axis[left : left + 3], values[left : left + 3], axis[i], width)
            b = _quadratic(axis[right : right + 3], values[right : right + 3], axis[i], width)
            coefficients.append([a[0], a[1] + b[0] - a[0], a[2] + b[1] - a[1], b[2] - a[2]])
        return coefficients
    # The reference falls back to its piecewise polynomial when fewer than
    # four supporting points are available for a cubic spline.
    if count < 4:
        return interval_coefficients(axis, values, 2)
    widths = [b - a for a, b in zip(axis, axis[1:])]
    rhs = [0.0] * count
    for i in range(1, count - 1):
        rhs[i] = 3 * (
            (values[i + 1] - values[i]) / widths[i] - (values[i] - values[i - 1]) / widths[i - 1]
        )
    diagonal, upper, solution = [0.0] * count, [0.0] * count, [0.0] * count
    if scheme == 4:
        rhs[0] = 3 * (values[1] - values[0]) / widths[0]
        rhs[-1] = -3 * (values[-1] - values[-2]) / widths[-1]
        diagonal[0], upper[0] = 2 * widths[0], 0.5
        solution[0] = rhs[0] / diagonal[0]
    else:
        diagonal[0] = 1.0
    for i in range(1, count - 1):
        diagonal[i] = 2 * (axis[i + 1] - axis[i - 1]) - widths[i - 1] * upper[i - 1]
        upper[i] = widths[i] / diagonal[i]
        solution[i] = (rhs[i] - widths[i - 1] * solution[i - 1]) / diagonal[i]
    second = [0.0] * count
    if scheme == 4:
        diagonal[-1] = widths[-1] * (2 - upper[-2])
        solution[-1] = (rhs[-1] - widths[-1] * solution[-2]) / diagonal[-1]
        second[-1] = solution[-1]
    coefficients = [None] * (count - 1)
    for i in range(count - 2, -1, -1):
        second[i] = solution[i] - upper[i] * second[i + 1]
        slope = (values[i + 1] - values[i]) / widths[i] - widths[i] * (
            second[i + 1] + 2 * second[i]
        ) / 3
        third = (second[i + 1] - second[i]) / (3 * widths[i])
        coefficients[i] = [
            values[i],
            slope * widths[i],
            second[i] * widths[i] ** 2,
            third * widths[i] ** 3,
        ]
    return coefficients


def parse_rpol_table(
    contents: str,
    *,
    default_dimension: int = 1,
    filename_options: str = "",
    api_options: str = "",
    default_boundary: str = "clamp",
) -> dict[str, Any]:
    """Read all three table formats; retain declared interpolation and boundary rules.

    API options override the header, and filename options override both.
    The reference scales logarithms after taking them.
    Column-count and allocation hints have no effect on the compiled kernel.
    """
    lines = contents.splitlines()
    header = next((line.strip() for line in lines if line.strip()), "")
    match = re.match(r"#@RPOL@(?:\[(.*?)\])?\s*([123])(?:\s|$)(.*)", header)
    dimension = int(match[2]) if match else default_dimension
    markup = match[1] if match else None
    text = " ".join((match[3] if match else "", api_options, filename_options))
    options = dict(boundary=default_boundary, scheme=1, x_log=False, y_log=False, value_log=False)
    scheme_given = False
    columns, scales = dict(x=0, y=1, z=2), dict(x=1.0, y=1.0, z=1.0)
    for token in re.split(r"[,\s]+", text.replace("OPTIONS:", "").strip()):
        if not token:
            continue
        key, _, value = token.lower().partition("=")
        if key == "clip":
            options["boundary"] = "clamp" if value in ("0", "off", "no") else "zero"
        elif key == "noclip":
            options["boundary"] = "clamp"
        elif key == "scheme":
            scheme_given = True
            try:
                options["scheme"] = int(value)
            except ValueError as error:
                raise TableImportError("invalid RPOL interpolation scheme") from error
            if not 0 <= options["scheme"] <= 4:
                raise TableImportError("unsupported RPOL interpolation scheme")
        elif key in ("xlog", "ylog", "zlog"):
            options[
                {
                    "xlog": "x_log",
                    "ylog": "y_log" if dimension != 1 else "value_log",
                    "zlog": "value_log",
                }[key]
            ] = True
        elif key in ("xcol", "ycol", "zcol"):
            try:
                columns[key[0]] = int(value) - 1
            except ValueError as error:
                raise TableImportError("invalid RPOL column index") from error
            if columns[key[0]] < 0:
                raise TableImportError("RPOL columns start at one")
            if dimension == 1 and key == "zcol":
                columns["y"] = columns["z"]
        elif key in ("xscale", "yscale", "zscale"):
            try:
                scales[key[0]] = {"deg2rad": math.pi / 180, "rad2deg": 180 / math.pi}.get(
                    value, float(value) if value not in ("deg2rad", "rad2deg") else 1
                )
            except ValueError as error:
                raise TableImportError("invalid RPOL scale") from error
            if not math.isfinite(scales[key[0]]):
                raise TableImportError("RPOL scale must be finite")
            if scales[key[0]] == 0:
                raise TableImportError("RPOL scale must not be zero")
        elif key in ("zxmax", "zxmin"):
            # Limits of the z range of the reference's plots: only the default, 0 (no limit), is
            # supported, as the table itself is not clipped by them.
            if value not in ("", "0"):
                raise TableImportError(f"unsupported RPOL option: {token}")
        elif key in ("verbose", "rows", "cols", "columns"):
            pass  # Output and allocation hints: no effect on the table.
        else:
            raise TableImportError(f"unknown RPOL option: {token}")
    rows = []
    for line in lines:
        line = line.strip()
        if markup and line.startswith(markup):
            line = line[len(markup) :].strip()
        elif not line or line.startswith(("#", "%")):
            continue
        try:
            row = [float(v) for v in line.split("#", 1)[0].split("%", 1)[0].split()]
        except ValueError as error:
            raise TableImportError("RPOL rows must be numeric") from error
        if not row or any(not math.isfinite(v) for v in row):
            raise TableImportError("RPOL rows must be finite")
        rows.append(row)

    def scaled(value, coordinate):
        logarithmic = options[{"x": "x_log", "y": "y_log", "z": "value_log"}[coordinate]]
        if logarithmic:
            if value <= 0:
                raise TableImportError("logarithmic RPOL values must be positive")
            return math.exp(math.log(value) * scales[coordinate])
        return value * scales[coordinate]

    try:
        if dimension == 1:
            samples = [
                (scaled(row[columns["x"]], "x"), scaled(row[columns["y"]], "z")) for row in rows
            ]
            axis = [x for x, _ in samples]
            values = [v for _, v in samples]
            if len(axis) < 2 or any(a >= b for a, b in zip(axis, axis[1:])):
                raise TableImportError("RPOL axis must increase strictly")
            if options["x_log"] and min(axis) <= 0 or options["value_log"] and min(values) <= 0:
                raise TableImportError("logarithmic RPOL values must be positive")
            x = [math.log(v) for v in axis] if options["x_log"] else axis
            y = [math.log(v) for v in values] if options["value_log"] else values
            options["coefficients"] = interval_coefficients(x, y, options["scheme"])
            return dict(x=axis, y=[], response=values, interpolation=options)
        if dimension == 2:
            if scheme_given and options["scheme"] != 1:
                raise TableImportError("RPOL 2-D tables are bilinear: only scheme 1 is supported")
            y = [scaled(v, "y") for v in rows[0]]
            if any(len(row) != len(y) + 1 for row in rows[1:]):
                raise TableImportError("RPOL matrix has invalid columns")
            samples = [
                (scaled(row[0], "x"), yy, scaled(value, "z"))
                for row in rows[1:]
                for yy, value in zip(y, row[1:])
            ]
        else:
            samples = [
                (
                    scaled(row[columns["x"]], "x"),
                    scaled(row[columns["y"]], "y"),
                    scaled(row[columns["z"]], "z"),
                )
                for row in rows
            ]
    except (IndexError, KeyError) as error:
        raise TableImportError("RPOL rows lack their requested columns") from error
    grid = {(x, y): value for x, y, value in samples}
    x, y = sorted({p[0] for p in samples}), sorted({p[1] for p in samples})
    if not grid or len(grid) != len(samples) or len(grid) != len(x) * len(y):
        raise TableImportError("RPOL table must be a complete grid with unique coordinates")
    if (
        options["x_log"]
        and min(x) <= 0
        or options["y_log"]
        and min(y) <= 0
        or options["value_log"]
        and min(grid.values()) <= 0
    ):
        raise TableImportError("logarithmic RPOL values must be positive")
    options["scheme"] = 1  # Reference implementation always uses bilinear 2-D interpolation.
    options["coefficients"] = []
    return dict(x=x, y=y, response=[grid[xx, yy] for xx in x for yy in y], interpolation=options)


def evaluate(
    axis: list[float], values: list[float], interpolation: dict[str, Any], x: float
) -> float:
    """Evaluate a 1-D RPOL table at x with its declared interpolation and boundary."""
    if x < axis[0] or x > axis[-1]:
        if interpolation["boundary"] == "zero":
            return 0.0
        x = min(max(x, axis[0]), axis[-1])
    transform = math.log if interpolation.get("x_log") else (lambda v: v)
    low = min(max(bisect.bisect_right(axis, x) - 1, 0), len(axis) - 2)
    t = (transform(x) - transform(axis[low])) / (transform(axis[low + 1]) - transform(axis[low]))
    scheme = interpolation.get("scheme", 1)
    logarithmic = interpolation.get("value_log", False)
    if scheme == 0:
        return values[low] if t < 0.5 else values[low + 1]
    if scheme >= 2:
        c = interpolation["coefficients"][low]
        value = ((c[3] * t + c[2]) * t + c[1]) * t + c[0]
    elif logarithmic:
        value = (1 - t) * math.log(values[low]) + t * math.log(values[low + 1])
    else:
        return (1 - t) * values[low] + t * values[low + 1]
    return math.exp(value) if logarithmic else value


_UNIT_TO_M = {"m": 1.0, "cm": 0.01, "mm": 0.001}
# sim_telarray mirror and camera-body shape codes.
SHAPES = {0: "circle", 1: "hexagon_flat_x", 2: "square", 3: "hexagon_flat_y"}


def number(value: object, name: str) -> float:
    """Return a finite number or raise."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ModelError(f"{name} must be a finite number")
    return float(value)


def length_unit_m(unit: object, name: str) -> float:
    """Metres per length `unit`."""
    if unit not in _UNIT_TO_M:
        raise ModelError(f"{name} has unsupported length unit {unit!r}")
    return _UNIT_TO_M[unit]


def signed_length_m(record: dict[str, Any], name: str) -> float:
    """Return a parameter record's length value in metres."""
    return number(record.get("value"), name) * length_unit_m(record.get("unit"), name)


def length_m(record: dict[str, Any], name: str, *, allow_zero: bool = False) -> float:
    """Return a positive (or non-negative) length in metres."""
    value = signed_length_m(record, name)
    if value < 0.0 or (value == 0.0 and not allow_zero):
        raise ModelError(f"{name} must be {'non-negative' if allow_zero else 'positive'}")
    return value


def components(record: dict[str, Any], name: str, count: int | None = None) -> tuple[list, list]:
    """Values and per-value units of a list parameter; a single unit applies to all values."""
    values, units = record.get("value"), record.get("unit")
    if not isinstance(values, list) or not values or count not in (None, len(values)):
        raise ModelError(f"{name} requires {count or 'a list of'} components")
    if units is None or isinstance(units, str):
        units = [units] * len(values)
    if not isinstance(units, list) or len(units) != len(values):
        raise ModelError(f"{name}: units do not match its values")
    return values, units


def lengths_m(record: dict[str, Any], name: str, count: int) -> list[float]:
    """Return a list of `count` lengths in metres."""
    values, units = components(record, name, count)
    return [number(v, name) * length_unit_m(u, name) for v, u in zip(values, units, strict=True)]


def _data_lines(contents: str, comments: tuple[str, ...] = ("#", "%")) -> list[list[str]]:
    rows = []
    for line in contents.splitlines():
        for marker in comments:
            line = line.split(marker, 1)[0]
        if line.strip():
            rows.append(line.split())
    return rows


def read_mirror_list(contents: str) -> list[dict[str, Any]]:
    """Read a sim_telarray mirror list: x, y, diameter, focal length [cm], shape, z [cm], and an
    optional panel id (a label of the mirror, a non-negative integer; the optics do not use it and
    the optical model does not store it).

    A focal length of zero asks for the automatic prescription; a negative one fixes its
    magnitude but asks for a random manufacturing error. A z of zero asks for the automatic
    height, a positive z fixes the height and a negative z gives its magnitude (see the compiler).
    """
    facets = []
    for fields in _data_lines(contents):
        if fields[:4] == ["mirror_x", "mirror_y", "mirror_diameter", "focal_length"]:
            continue  # ECSV column header
        if len(fields) not in (6, 7):
            raise ModelError(
                f"mirror list row {len(facets)}: expected 6 or 7 columns (x, y, diameter, "
                "focal length, shape, z and an optional panel id)"
            )
        try:
            x, y, diameter, focal, shape, z = (float(v) for v in fields[:6])
            panel = float(fields[6]) if len(fields) == 7 else None
        except ValueError as error:
            raise ModelError(f"mirror list row {len(facets)}: invalid number") from error
        if not all(math.isfinite(v) for v in (x, y, diameter, focal, shape, z)):
            raise ModelError(f"mirror list row {len(facets)}: non-finite value")
        if panel is not None and (not math.isfinite(panel) or panel < 0 or panel != int(panel)):
            raise ModelError(f"mirror list row {len(facets)}: the panel id must be an integer >= 0")
        if shape not in SHAPES or diameter <= 0.0:
            raise ModelError(f"mirror list row {len(facets)}: invalid shape or diameter")
        facets.append({
            "id": len(facets),
            "centre_m": [0.01 * x, 0.01 * y, 0.01 * z],
            "shape": SHAPES[int(shape)],
            "diameter_m": 0.01 * diameter,
            "catalogue_focal_length_m": 0.01 * focal,
            "panel_id": None if panel is None else int(panel),
        })
    if not facets:
        raise ModelError("mirror list contains no mirrors")
    return facets


def read_segmentation(contents: str) -> list[dict[str, Any]]:
    """Read sim_telarray segmentation records (hex, yhex, square, circle, ring, polygon)."""
    names = {
        "hexsegments": "hex",
        "yhexsegments": "yhex",
        "ringsegments": "ring",
        "squaresegments": "square",
        "circularsegments": "circle",
        "polygonsegments": "polygon",
    }
    segments: list[dict[str, Any]] = []
    for row, fields in enumerate(_data_lines(contents, ("#",)), start=1):
        token = fields[0].lower()
        kind = next((value for name, value in names.items() if name.startswith(token)), token)
        try:
            count = int(fields[1])
            values = [float(v) for v in " ".join(fields[2:]).replace(",", " ").split()]
        except (ValueError, IndexError) as error:
            raise ModelError(f"segmentation row {row}: invalid number") from error
        if count < 1 or not all(math.isfinite(v) for v in values):
            raise ModelError(f"segmentation row {row}: invalid count or value")
        if kind == "polygon":
            if count != 1 or len(values) < 7 or len(values) % 2 == 0:
                raise ModelError(f"segmentation row {row}: polygon needs an angle and points")
            angle = math.radians(values[0])
            vertices = [
                [
                    0.01 * (x * math.cos(angle) + y * math.sin(angle)),
                    0.01 * (-x * math.sin(angle) + y * math.cos(angle)),
                ]
                for x, y in zip(values[1::2], values[2::2])
            ]
            if vertices[0] == vertices[-1]:
                vertices.pop()
            turns = [
                (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
                for a, b, c in zip(
                    vertices, vertices[1:] + vertices[:1], vertices[2:] + vertices[:2]
                )
            ]
            if len(vertices) < 3 or not (all(t > 0 for t in turns) or all(t < 0 for t in turns)):
                raise ModelError(f"segmentation row {row}: polygon must be strictly convex")
            segments.append({"id": len(segments), "shape": "polygon", "vertices_xy_m": vertices})
        elif kind in ("hex", "yhex", "square", "circle"):
            if count != 1 or len(values) not in (3, 4) or values[2] <= 0:
                raise ModelError(f"segmentation row {row}: invalid footprint")
            x, y, diameter, rotation = (*values, 0.0)[:4]
            segments.append({
                "id": len(segments),
                "shape": {"square": "square", "circle": "circle"}.get(kind, "hexagon"),
                "centre_xy_m": [0.01 * x, 0.01 * y],
                "diameter_m": 0.01 * diameter,
                "rotation_deg": rotation + (90.0 if kind == "yhex" else 0.0),
            })
        elif kind == "ring":
            if len(values) not in (3, 4, 5):
                raise ModelError(f"segmentation row {row}: invalid ring")
            inner, outer, span, start, gap = (*values, 0.0, 0.0)[:5]
            if inner < 0 or outer <= inner or span <= 0 or gap < 0 or count * span > 360 + 1e-9:
                raise ModelError(f"segmentation row {row}: invalid ring")
            first = len(segments)
            segments.extend(
                {
                    "id": first + i,
                    "shape": "annular_sector",
                    "inner_radius_m": 0.01 * inner,
                    "outer_radius_m": 0.01 * outer,
                    "start_deg": start + i * span,
                    "span_deg": span,
                    "gap_m": 0.01 * gap,
                    "gap_at_start": False,
                }
                for i in range(count)
            )
        else:
            raise ModelError(f"segmentation row {row}: unsupported type {fields[0]}")
    if not segments:
        raise ModelError("segmentation contains no segments")
    return segments


def segmentation_from_record(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Convert a structured simulation-models segmentation into sim_telarray rows."""
    groups = record.get("value")
    if not isinstance(groups, list) or not groups:
        raise ModelError("structured segmentation must list footprint groups")
    lines = []
    for group in groups:
        kind = group.get("kind") if isinstance(group, dict) else None
        footprint = kind in ("hex", "yhex", "square", "circle")
        names = (
            ("x", "y", "diameter", "rotation")
            if footprint
            else ("r_min", "r_max", "dphi", "phi0", "gap")
        )
        optional = {"rotation"} if footprint else {"phi0", "gap"}
        if kind not in ("hex", "yhex", "square", "circle", "ring") or set(group) - {
            "kind",
            "count",
            *names,
        }:
            raise ModelError("invalid structured segmentation fields")

        def quantity(name: str, group: dict = group, optional: set = optional) -> float:
            value = group.get(name)
            if value is None and name in optional:
                return 0.0
            if not isinstance(value, dict) or set(value) != {"value", "unit"}:
                raise ModelError(f"segmentation {name} requires value and unit")
            if name in ("rotation", "dphi", "phi0"):
                if value["unit"] != "deg":
                    raise ModelError(f"segmentation {name} requires degrees")
                return number(value["value"], name)
            return signed_length_m(value, name) * 100

        if isinstance(group.get("count"), bool) or not isinstance(group.get("count"), int):
            raise ModelError("segmentation count must be an integer")
        values = " ".join(format(quantity(name), ".17g") for name in names)
        lines.append(f"{kind} {group['count']} {values}")
    return read_segmentation("\n".join(lines))


def secondary_segment_frame(segments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Rotate the secondary mirror frame by pi about y."""
    result = []
    for segment in segments:
        if segment["shape"] == "polygon":
            result.append({
                **segment,
                "vertices_xy_m": [[-x, y] for x, y in segment["vertices_xy_m"]],
            })
        elif segment["shape"] == "annular_sector":
            start = 180 - segment["start_deg"] - segment["span_deg"]
            result.append({**segment, "start_deg": start, "gap_at_start": True})
        else:
            x, y = segment["centre_xy_m"]
            result.append({
                **segment,
                "centre_xy_m": [-x, y],
                "rotation_deg": 180 - segment["rotation_deg"],
            })
    return result


def read_response(
    contents: str, name: str, value_column: str = "reflectivity", *, filename_options: str = ""
) -> dict[str, Any]:
    """Read a response versus wavelength, optionally versus incidence angle.

    Returns ``{"x": wavelengths_nm, "y": incidence_deg (empty for 1-D),
    "values": [x][y] row-major, "interpolation": RPOL options}``.
    """
    rows = _data_lines(contents)
    if not rows:
        raise ModelError(f"{name} is empty")
    try:
        float(rows[0][0])
        numeric_first = True
    except ValueError, IndexError:
        numeric_first = False
    if contents.lstrip().startswith("#@RPOL@") or numeric_first:
        try:
            table = parse_rpol_table(
                contents, api_options="clip=0,yscale=deg2rad", filename_options=filename_options
            )
        except TableImportError as error:
            raise ModelError(f"{name}: {error}") from error
        result = {
            "x": table["x"],
            "y": [math.degrees(v) for v in table["y"]],
            "values": table["response"],
            "interpolation": table["interpolation"],
        }
    else:
        header = rows[0]
        columns = {"wavelength", value_column} | (
            {"incidence_angle"} if "incidence_angle" in header else set()
        )
        if not columns.issubset(header) or len(header) != len(set(header)):
            raise ModelError(f"{name} lacks wavelength and {value_column} columns")
        samples: dict[tuple[float, float], list[float]] = {}
        for row in rows[1:]:
            if len(row) != len(header):
                raise ModelError(f"{name}: wrong column count")
            entry = dict(zip(header, row, strict=True))
            try:
                key = (float(entry["wavelength"]), float(entry.get("incidence_angle", 0.0)))
                samples.setdefault(key, []).append(float(entry[value_column]))
            except ValueError as error:
                raise ModelError(f"{name}: invalid number") from error
        x = sorted({k[0] for k in samples})
        y = sorted({k[1] for k in samples}) if "incidence_angle" in header else []
        if len(samples) != len(x) * max(1, len(y)):
            raise ModelError(f"{name} requires a complete wavelength-angle grid")
        # Repeated measurements of one wavelength are averaged.
        values = [sum(samples[xx, yy]) / len(samples[xx, yy]) for xx in x for yy in (y or [0.0])]
        result = {
            "x": x,
            "y": y,
            "values": values,
            "interpolation": {
                "boundary": "clamp",
                "scheme": 1,
                "x_log": False,
                "value_log": False,
                "coefficients": [],
            },
        }
    if len(result["x"]) < 2 or any(not 0.0 <= v <= 1.0 for v in result["values"]):
        raise ModelError(f"{name} needs at least two wavelengths and values in [0, 1]")
    return result


def read_spatial_map(contents: str, name: str, *, filename_options: str = "") -> dict[str, Any]:
    """Read a 2-D relative response map in cm coordinates; return metres."""
    try:
        table = parse_rpol_table(contents, default_dimension=3, filename_options=filename_options)
    except TableImportError as error:
        raise ModelError(f"{name}: {error}") from error
    if not table["y"] or any(not 0.0 <= v <= 1.0 for v in table["response"]):
        raise ModelError(f"{name} must be a 2-D table of fractions")
    ny = len(table["y"])
    return {
        "x_m": [0.01 * v for v in table["x"]],
        "y_m": [0.01 * v for v in table["y"]],
        "values": [table["response"][i * ny : (i + 1) * ny] for i in range(len(table["x"]))],
        "outside": "zero" if table["interpolation"]["boundary"] == "zero" else "clamp",
    }


_STRUCTURE_TEXT = ("id", "group", "shape")


def read_structure_table(contents: str, name: str) -> dict[str, Any]:
    """Read an ECSV structure table; check units (m) and the coordinate frame.

    Returns ``{"frame": str, "columns": [str], "rows": [{column: value}]}``;
    numeric columns are floats, ``nan`` marking a column the row does not use.
    """
    header = None
    units: dict[str, str] = {}
    frame = None
    rows = []
    for line in contents.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            column = re.match(
                r"#\s*-\s*\{name:\s*(\w+),\s*(?:unit:\s*(\w+),\s*)?datatype", stripped
            )
            if column and column[2]:
                units[column[1]] = column[2]
            if match := re.match(r"#\s*coordinate_frame:\s*(\S+)", stripped):
                frame = match[1]
            continue
        if not stripped:
            continue
        fields = stripped.split()
        if header is None:
            header = fields
            continue
        if len(fields) != len(header):
            raise ModelError(f"{name} row {len(rows)}: wrong column count")
        row: dict[str, Any] = {}
        for column, field in zip(header, fields, strict=True):
            if column in _STRUCTURE_TEXT:
                row[column] = field
                continue
            try:
                row[column] = float(field)
            except ValueError as error:
                raise ModelError(f"{name} row {len(rows)}: invalid {column}") from error
            if math.isinf(row[column]):
                raise ModelError(f"{name} row {len(rows)}: infinite {column}")
        rows.append(row)
    if header is None or not rows:
        raise ModelError(f"{name} is empty")
    for column in header:
        if column not in _STRUCTURE_TEXT and units.get(column) != "m":
            raise ModelError(f"{name}: column {column} must be declared in metres")
    if frame is None:
        raise ModelError(f"{name} does not declare a coordinate_frame")
    return {"frame": frame, "columns": header, "rows": rows}
