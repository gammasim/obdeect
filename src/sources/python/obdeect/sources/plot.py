"""Plots of the emitted profiles of a light source (optional: needs ``dlamb[plot]``).

The plots sample photons with ``dlamb.source.emit``: they show what a source emits,
not the analytic profile. matplotlib is imported when a function is called.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from obdeect.sources.source import CHUNK_ROWS, build_source, emit_chunk


def _sample(spec: Mapping[str, Any], seed: int, photons: int) -> tuple[dict, np.ndarray]:
    """Up to `photons` emitted photons and their weights; the first chunk only."""
    count = build_source(spec).photon_count(seed, 0, 0)
    rows = min(photons, count, CHUNK_ROWS)
    if rows == 0:
        raise ValueError("the source emits no photons")
    chunk = emit_chunk(spec, seed=seed, rows=rows)
    return chunk, chunk["w"]


def _angle_from_axis_deg(chunk: dict, spec: Mapping[str, Any]) -> np.ndarray:
    """Angle [deg] between each emitted direction and the emitter axis (or the star's source)."""
    if "emitter" in spec:
        axis = np.asarray(spec["emitter"]["pointing"], dtype=float)
        axis = axis / np.linalg.norm(axis)
    else:
        towards = np.asarray(spec["star"]["towards_source"], dtype=float)
        axis = -towards / np.linalg.norm(towards)
    cosine = chunk["cx"] * axis[0] + chunk["cy"] * axis[1] + chunk["cz"] * axis[2]
    return np.degrees(np.arccos(np.clip(cosine.astype(float), -1.0, 1.0)))


def profile_figure(spec: Mapping[str, Any], *, seed: int = 1, photons: int = 200_000):
    """Figure with three panels: angle to the axis, wavelength and emission time.

    Histograms are weighted by the photon weights (disc targets). Returns the
    matplotlib figure; the caller saves or shows it.
    """
    import matplotlib.pyplot as plt

    chunk, weights = _sample(spec, seed, photons)
    figure, (angle_axis, wavelength_axis, time_axis) = plt.subplots(1, 3, figsize=(12, 3.6))
    angle_axis.hist(_angle_from_axis_deg(chunk, spec), bins=60, weights=weights, color="#4C78A8")
    angle_axis.set_xlabel("angle to axis [deg]")
    wavelength_axis.hist(chunk["wavelength"], bins=60, weights=weights, color="#F58518")
    wavelength_axis.set_xlabel("wavelength [nm]")
    time_axis.hist(chunk["t"], bins=60, weights=weights, color="#54A24B")
    time_axis.set_xlabel("emission time [ns]")
    for axis in (angle_axis, wavelength_axis, time_axis):
        axis.set_ylabel("photons (weighted)")
    figure.tight_layout()
    return figure
