"""Atmospheric profile, refraction and light timing on numpy arrays.

Heights are in m a.s.l., positions in m, times in ns. Directions are the direction cosines
(cx, cy) of a downward photon (cz = sqrt(1 - cx^2 - cy^2) > 0). See docs/refraction.md.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from dlufd import _core


def _as_float_array(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float64)


class AtmosphereProfile:
    """CORSIKA/sim_telarray atmprof table with the atmo.c interpolation of density,
    thickness and n-1 (see docs/refraction.md)."""

    def __init__(self, path: str | Path):
        self._core = _core.AtmosphereProfile(str(path))

    @property
    def levels(self) -> int:
        return self._core.levels

    @property
    def bottom_m(self) -> float:
        return self._core.bottom_m

    @property
    def top_m(self) -> float:
        return self._core.top_m

    def density(self, height_m):
        """Density [g/cm^3] at heights [m a.s.l.]; same shape as the input."""
        return self._apply(self._core.density, height_m)

    def thickness(self, height_m):
        """Atmospheric thickness [g/cm^2] at heights [m a.s.l.]; same shape as the input."""
        return self._apply(self._core.thickness, height_m)

    def refractive_index_minus_one(self, height_m):
        """Refractive index minus one at heights [m a.s.l.]; same shape as the input."""
        return self._apply(self._core.refractive_index_minus_one, height_m)

    def group_index(self, observation_level_m: float) -> float:
        """Refractive index at the observation level (obdeect ambient_group_index)."""
        return self._core.group_index(float(observation_level_m))

    @staticmethod
    def _apply(function, values):
        array = _as_float_array(values)
        result = function(np.ascontiguousarray(array.reshape(-1)))
        return result.reshape(array.shape)


class RefractionTables:
    """Refraction tables of CORSIKA atmo.c for a profile and an observation level."""

    def __init__(self, profile: AtmosphereProfile, observation_level_m: float):
        self._core = _core.RefractionTables(profile._core, float(observation_level_m))

    @property
    def etadsn(self) -> float:
        return self._core.etadsn

    @property
    def observation_level_m(self) -> float:
        return self._core.observation_level_m


def refract(tables: RefractionTables, height_m, cx, cy, x_m, y_m, t_ns):
    """Bend photons emitted at heights [m a.s.l.] to the observation level.

    Inputs are one-dimensional arrays of equal length: emission height, direction cosines
    (cx, cy), straight-line position (x_m, y_m) and time (t_ns) at the observation level.
    Returns (cx, cy, x_m, y_m, t_ns, bent) where bent is a bool array; photons with bent False
    (direction outside the unit circle) are returned unchanged.
    """
    arrays = [
        np.ascontiguousarray(_as_float_array(values))
        for values in (height_m, cx, cy, x_m, y_m, t_ns)
    ]
    if any(array.ndim != 1 for array in arrays):
        raise ValueError("refract: all inputs must be one-dimensional")
    if len({array.shape[0] for array in arrays}) != 1:
        raise ValueError("refract: all inputs must have the same length")
    new_cx, new_cy, new_x, new_y, new_t, bent = _core.refract(tables._core, *arrays)
    return new_cx, new_cy, new_x, new_y, new_t, bent.astype(bool)


def air_delay_ns(
    profile: AtmosphereProfile, height_low_m: float, height_high_m: float, zenith_rad: float
) -> float:
    """Extra light travel time along a straight path between two heights [ns]."""
    return _core.straight_path_air_delay_ns(
        profile._core, float(height_low_m), float(height_high_m), float(zenith_rad)
    )


def propagation_time_ns(path_m: float, group_index: float) -> float:
    """Light travel time over a path length in a medium of the given index [ns]."""
    return _core.propagation_time_ns(float(path_m), float(group_index))
