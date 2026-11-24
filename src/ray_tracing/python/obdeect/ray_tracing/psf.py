"""Point-spread function and effective area from traced photons."""

from __future__ import annotations

from typing import Any

import numpy as np
from obdeect.ray_tracing.model import DETECTED, STATUS_NAMES


def containment(
    x: np.ndarray, y: np.ndarray, weights: np.ndarray | None = None, fraction: float = 0.8
) -> dict:
    """Centroid and the radius around it containing `fraction` of the weight.

    Raises ValueError for a fraction outside (0, 1], for no photons, and for weights that sum
    to zero: the radius is then undefined, and no other weighting is substituted.
    """
    if not 0.0 < fraction <= 1.0:
        raise ValueError("fraction must lie in (0, 1]")
    if len(x) == 0:
        raise ValueError("no photons on the focal surface")
    w = np.ones_like(x) if weights is None else np.asarray(weights, dtype=float)
    total = float(w.sum())
    if not total > 0.0:
        raise ValueError("the weights on the focal surface sum to zero: the PSF is undefined")
    cx, cy = np.dot(w, x) / total, np.dot(w, y) / total
    r = np.hypot(x - cx, y - cy)
    order = np.argsort(r)
    cumulative = np.cumsum(w[order])
    index = min(int(np.searchsorted(cumulative, fraction * total)), len(r) - 1)
    return {
        "centroid_x_m": float(cx),
        "centroid_y_m": float(cy),
        "radius_m": float(r[order][index]),
    }


def summary(
    arrivals: dict[str, np.ndarray],
    sampling_area_m2: float,
    fraction: float = 0.8,
    *,
    sampled: int | None = None,
) -> dict[str, Any]:
    """Status counts, effective area and focal-plane PSF of one traced source.

    The effective area is the sum of the optical weights over the sampled photons, times the
    sampled area. `sampled` is the number of photons sampled: it is required for results of
    ``detected_only`` tracing (identified by the ``index`` column), whose arrays hold the
    detected photons only. Without detected photons, or with weights that sum to zero, no PSF
    is given; ``detected_weight`` reports the weight on the focal surface.
    """
    status = arrivals["status"]
    if sampled is None:
        if "index" in arrivals:
            raise ValueError("detected-only results need `sampled`, the number of photons sampled")
        sampled = len(status)
    if sampled < len(status) or sampled < 1:
        raise ValueError("sampled must be at least 1 and at least the number of arrivals")
    detected = status == DETECTED
    counts = {name: int(np.count_nonzero(status == i)) for i, name in enumerate(STATUS_NAMES)}
    result: dict[str, Any] = {
        "photons": int(sampled),
        "status": {k: v for k, v in counts.items() if v},
        "effective_area_m2": float(arrivals["weight"].sum() / sampled * sampling_area_m2),
        "detected_weight": float(arrivals["weight"][detected].sum()),
    }
    if detected.any() and result["detected_weight"] > 0.0:
        position = arrivals["position_m"][detected]
        psf = containment(position[:, 0], position[:, 1], arrivals["weight"][detected], fraction)
        result.update(psf)
        result[f"d{round(100 * fraction)}_m"] = 2 * psf.pop("radius_m")
    return result
