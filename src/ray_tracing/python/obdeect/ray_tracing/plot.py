"""Diagnostic plot: the telescope seen from the sky, and the focal-plane image.

Matplotlib is optional (the ``plot`` extra). The figure is drawn on its own canvas, so importing
this module does not select a matplotlib backend for the caller.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from obdeect.ray_tracing.model import DETECTED, OpticalModel

_CORNERS = {"circle": 48, "square": 4, "hexagon_flat_x": 6, "hexagon_flat_y": 6}


def _outline(shape: str, half_width: float) -> np.ndarray:
    """Aperture outline in its (u, v) frame."""
    n = _CORNERS[shape]
    if shape == "circle":
        angles, radius = np.linspace(0, 2 * np.pi, n + 1), half_width
    elif shape == "square":
        angles, radius = np.pi / 4 + np.arange(n + 1) * np.pi / 2, half_width * math.sqrt(2)
    else:  # Corners on the axis parallel to the flats.
        start = np.pi / 2 if shape == "hexagon_flat_x" else 0.0
        angles, radius = start + np.arange(n + 1) * np.pi / 3, half_width * 2 / math.sqrt(3)
    return radius * np.column_stack([np.cos(angles), np.sin(angles)])


def _in_plane(outline: np.ndarray, centre: np.ndarray, normal: np.ndarray, tangent: np.ndarray):
    """Outline points of an aperture in its own (u, v) frame, as 3-D points."""
    u = np.array(tangent, dtype=float)
    u /= np.linalg.norm(u)
    v = np.cross(np.array(normal, dtype=float), u)
    v /= np.linalg.norm(v)
    return np.asarray(centre, dtype=float) + outline[:, :1] * u + outline[:, 1:] * v


def plot(model: OpticalModel, arrivals: dict[str, np.ndarray] | None, output: Path) -> None:
    optics = model.document["optics"]
    figure = Figure(figsize=(12 if arrivals is not None else 6, 6))
    FigureCanvasAgg(figure)
    sky = figure.add_subplot(1, 2 if arrivals is not None else 1, 1)
    primary = optics["primary"]
    for facet in primary.get("facets", []):
        points = _in_plane(
            _outline(facet["shape"], 0.5 * facet["diameter_m"]),
            facet["centre_m"],
            facet["normal"],
            facet["tangent"],
        )
        sky.fill(points[:, 0], points[:, 1], color="0.85", ec="0.4", lw=0.5)
    if "asphere" in primary:
        for radius in (primary["asphere"]["inner_radius_m"], primary["asphere"]["outer_radius_m"]):
            sky.add_patch(Circle((0, 0), radius, fill=False, color="0.4"))
    for obscurer in optics["obscurers"]:
        if obscurer["shape"] == "cylinder":
            ends = np.array([obscurer["start_m"], obscurer["end_m"]])
            sky.plot(ends[:, 0], ends[:, 1], color="C3", lw=1)
        elif obscurer["shape"] == "frustum":
            # Both ends as circles (the outer radius, wall included), seen along the axis.
            for point, radius in zip(
                (obscurer["start_m"], obscurer["end_m"]),
                (obscurer["start_radius_m"], obscurer["end_radius_m"]),
                strict=True,
            ):
                sky.add_patch(
                    Circle(
                        (point[0], point[1]),
                        radius + obscurer["thickness_m"],
                        fill=False,
                        color="C3",
                        lw=1,
                    )
                )
        elif obscurer["shape"] == "polyhedron":
            vertices = np.array(obscurer["vertices_m"])
            sky.plot(vertices[:, 0], vertices[:, 1], ".", color="C3", ms=2)
        elif obscurer["shape"] == "plate":
            outline = _in_plane(
                _outline(obscurer["aperture"], 0.5 * obscurer["diameter_m"]),
                obscurer["centre_m"],
                obscurer["normal"],
                obscurer["tangent"],
            )
            sky.plot(outline[:, 0], outline[:, 1], color="C0", lw=1)
    sky.set(aspect="equal", xlabel="x [m]", ylabel="y [m]", title="Telescope seen from the sky")
    if arrivals is not None:
        detected = arrivals["status"] == DETECTED
        position = arrivals["position_m"][detected]
        image = figure.add_subplot(1, 2, 2)
        image.hist2d(
            1e3 * position[:, 0],
            1e3 * position[:, 1],
            bins=100,
            weights=arrivals["weight"][detected],
            cmap="viridis",
        )
        image.set(aspect="equal", xlabel="x [mm]", ylabel="y [mm]", title="Focal-plane image")
    figure.suptitle(f"optical model {model.sha256[:12]}")
    figure.tight_layout()
    figure.savefig(output)
