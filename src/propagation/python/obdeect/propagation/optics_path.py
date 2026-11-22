"""In-telescope path length of an optical model.

The extinction inside a telescope covers the path from the primary reflector vertex to the
focus. The optical model gives the positions on the optical axis.
"""

from __future__ import annotations

import argparse
import os
from collections.abc import Sequence
from pathlib import Path


def in_telescope_path_m(optical_model) -> float:
    """Return the in-telescope path length [m] of an optical model.

    Definition, with all positions on the optical axis in the telescope frame [m]:

    - single mirror: ``|f0 - p0|``, the focal length from the primary vertex p0 to the focus f0;
    - dual mirror: ``|p0 - s0| + |s0 - f0|``, from the primary vertex p0 via the secondary vertex
      s0 to the focus f0.

    The primary vertex is the dish vertex at ``-mirror_offset`` (for a segmented primary the
    compiler stores it as ``optics.primary.vertex_z_m``). The focus f0 is the focal surface on the
    axis (plane height, or the vertex of an asphere).

    Parameters
    ----------
    optical_model
        A loaded ``OpticalModel`` or the path of an optical model file (JSON).

    Raises
    ------
    ValueError
        If the model stores no primary vertex (models compiled before ``vertex_z_m`` was written
        for segmented primaries: compile the model again).
    TypeError
        If the argument is neither an optical model nor a path.
    """
    from obdeect.ray_tracing.model import OpticalModel

    if isinstance(optical_model, OpticalModel):
        model = optical_model
    elif isinstance(optical_model, (str, os.PathLike)):
        model = OpticalModel.load(Path(optical_model))
    else:
        raise TypeError("optical_model must be an OpticalModel or a path to an optical model file")

    optics = model.document["optics"]
    primary_z = _vertex_z_m(optics["primary"], "primary")
    secondary = optics.get("secondary")
    secondary_z = None if secondary is None else _vertex_z_m(secondary, "secondary")
    return _path_m(primary_z, secondary_z, model.focal_surface_z_m)


def _vertex_z_m(mirror: dict, name: str) -> float:
    """Vertex height [m] of a mirror: the asphere vertex, or the stored vertex of facets."""
    if "vertex_z_m" in mirror:
        return float(mirror["vertex_z_m"])
    if "asphere" in mirror:
        return float(mirror["asphere"]["vertex_z_m"])
    raise ValueError(
        f"optical model has no {name} vertex: a segmented {name} without vertex_z_m was compiled "
        "before the dish vertex was stored; compile the model again"
    )


def _path_m(primary_z: float, secondary_z: float | None, focus_z: float) -> float:
    """Path [m] from the primary vertex to the focus, via the secondary vertex if present."""
    if secondary_z is None:
        return abs(focus_z - primary_z)
    return abs(primary_z - secondary_z) + abs(secondary_z - focus_z)


def main(argv: Sequence[str] | None = None) -> int:
    """Command line: ``obdeect propagation pathlen MODEL``. Prints the path in metres."""
    parser = argparse.ArgumentParser(
        prog="obdeect propagation pathlen",
        description="Print the in-telescope path length [m] of an optical model (JSON).",
    )
    parser.add_argument("model", help="optical model file (JSON)")
    args = parser.parse_args(argv)
    try:
        path = in_telescope_path_m(args.model)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"obdeect propagation pathlen: {error}\n")
    print(f"{path:.6f}")
    return 0
