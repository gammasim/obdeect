"""Conservative bounds of optical elements, from the optics records of an optical model.

A bound is an upper bound (of a height or of a dot product with a direction) over the whole
element, not only over its vertices: asphere sags, facet footprints and the walls of obscurers
included. The compiler places the star-light disk with `optics_top_z_m`; the model checks input
planes with the same bounds, so both agree on where the optics end.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np
from obdeect.ray_tracing.errors import ModelError

# Radial samples of an asphere mirror's aperture for its support.
ASPHERE_SAMPLES = 2049
UP = np.array([0.0, 0.0, 1.0])


def circumradius(aperture: str, diameter_m: float) -> float:
    """Radius [m] of the circle that encloses an aperture of the given diameter (flat to flat
    for hexagons, corner to corner for squares)."""
    half_width = 0.5 * float(diameter_m)
    if aperture == "circle":
        return half_width
    if aperture == "square":
        return half_width * math.sqrt(2.0)
    if aperture in ("hexagon_flat_x", "hexagon_flat_y"):
        return half_width * 2.0 / math.sqrt(3.0)
    raise ModelError(f"unknown aperture shape {aperture!r}")


def end_radii(obscurer: Mapping[str, Any]) -> tuple[float, float]:
    """Outer radii [m] at the start and end of a cylinder or frustum (wall included)."""
    if "start_radius_m" in obscurer:
        thickness = float(obscurer["thickness_m"])
        return (
            float(obscurer["start_radius_m"]) + thickness,
            float(obscurer["end_radius_m"]) + thickness,
        )
    return (0.5 * float(obscurer["diameter_m"]),) * 2


def obscurer_top_z_m(obscurer: Mapping[str, Any]) -> float:
    """Highest z [m] of one obscurer record."""
    if "vertices_m" in obscurer:
        return max(float(v[2]) for v in obscurer["vertices_m"])
    if "centre_m" in obscurer:
        normal = np.asarray(obscurer["normal"], dtype=float)
        normal /= np.linalg.norm(normal)
        radius = circumradius(obscurer["aperture"], obscurer["diameter_m"])
        return float(obscurer["centre_m"][2]) + radius * math.sqrt(max(0.0, 1.0 - normal[2] ** 2))
    start, end = (np.asarray(obscurer[key], dtype=float) for key in ("start_m", "end_m"))
    axis = (end - start) / np.linalg.norm(end - start)
    sine = math.sqrt(max(0.0, 1.0 - axis[2] ** 2))
    return max(
        float(point[2]) + radius * sine
        for point, radius in zip((start, end), end_radii(obscurer), strict=True)
    )


def obscurer_support(obscurer: Mapping[str, Any], direction: np.ndarray) -> float:
    """Maximum dot product with a direction over an obscurer's geometry."""
    if "vertices_m" in obscurer:
        return max(float(np.dot(v, direction)) for v in obscurer["vertices_m"])
    if "centre_m" in obscurer:
        normal = np.asarray(obscurer["normal"], dtype=float)
        normal /= np.linalg.norm(normal)
        radius = circumradius(obscurer["aperture"], obscurer["diameter_m"])
        radial = math.sqrt(
            max(0.0, float(np.dot(direction, direction) - np.dot(direction, normal) ** 2))
        )
        return float(np.dot(obscurer["centre_m"], direction)) + radius * radial
    start, end = (np.asarray(obscurer[key], dtype=float) for key in ("start_m", "end_m"))
    axis = (end - start) / np.linalg.norm(end - start)
    radial = math.sqrt(max(0.0, float(np.dot(direction, direction) - np.dot(direction, axis) ** 2)))
    return max(
        float(np.dot(point, direction)) + radius * radial
        for point, radius in zip((start, end), end_radii(obscurer), strict=True)
    )


def asphere_support(
    asphere: Mapping[str, Any], direction: np.ndarray, *, margin: bool = True
) -> float:
    """Bound on the dot product over an asphere (a mirror or the focal surface).

    An asphere point at radius r has the direction's dot product hypot(d_x, d_y) r + d_z z(r).
    The maximum over the samples is taken; with `margin`, the interpolation error between radial
    samples is added, bounded with the maximum absolute second derivative of that polynomial.
    """
    radii = np.linspace(
        float(asphere["inner_radius_m"]), float(asphere["outer_radius_m"]), ASPHERE_SAMPLES
    )
    scale = float(asphere["radial_scale_m"])
    z = np.full(radii.shape, float(asphere["vertex_z_m"]))
    curvature_bound = 0.0
    magnitude = abs(float(asphere["vertex_z_m"]))
    outer_scaled = float(asphere["outer_radius_m"]) / scale
    for i, coefficient in enumerate(asphere["coefficients_m"]):
        coefficient = float(coefficient)
        z += coefficient * (radii / scale) ** (2 * i)
        magnitude += abs(coefficient) * outer_scaled ** (2 * i)
        if i:
            curvature_bound += (
                abs(coefficient) * (2 * i) * (2 * i - 1) * outer_scaled ** (2 * i - 2) / scale**2
            )
    radial = math.hypot(float(direction[0]), float(direction[1]))
    dz = abs(float(direction[2]))
    step = (float(asphere["outer_radius_m"]) - float(asphere["inner_radius_m"])) / (
        ASPHERE_SAMPLES - 1
    )
    rounding = (
        np.finfo(float).eps
        * (16 + 4 * len(asphere["coefficients_m"]))
        * (radial * float(asphere["outer_radius_m"]) + dz * magnitude)
    )
    peak = float(np.max(radial * radii + float(direction[2]) * z))
    if not margin:
        return peak
    return peak + dz * curvature_bound * step**2 / 8 + rounding


def mirror_support(
    mirror: Mapping[str, Any], direction: np.ndarray, *, margin: bool = True
) -> float:
    """Conservative bound on the dot product over a mirror's surface.

    A facet's footprint lies within its circumradius plus its sag about the centre, so
    that bound is used. `margin` is passed to `asphere_support`.
    """
    if "facets" in mirror:
        best = -math.inf
        for facet in mirror["facets"]:
            radius = float(facet.get("radius_of_curvature_m", 0.0))
            reach = circumradius(facet["shape"], facet["diameter_m"])
            sag = radius - math.sqrt(max(0.0, radius * radius - reach * reach)) if radius else 0.0
            best = max(best, float(np.dot(facet["centre_m"], direction)) + reach + sag)
        return best
    return asphere_support(mirror["asphere"], direction, margin=margin)


def optics_top_z_m(optics: Mapping[str, Any], *, margin: bool = True) -> float:
    """Highest z [m] of any optical element of an optics record: the mirrors (primary and
    secondary), the focal surface and every obscurer on any leg. `margin` as in
    `asphere_support` (the interpolation margin of the asphere samples)."""
    tops = [
        mirror_support(optics[role], UP, margin=margin)
        for role in ("primary", "secondary")
        if optics.get(role) is not None
    ]
    focal = optics["focal_surface"]
    tops.append(
        asphere_support(focal["asphere"], UP, margin=margin)
        if "asphere" in focal
        else float(focal["plane_z_m"])
    )
    tops += [obscurer_top_z_m(o) for o in optics["obscurers"]]
    return max(tops)
