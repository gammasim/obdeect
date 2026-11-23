"""Direction and value helpers of the presets (laser-beam, lidar, night-sky).

Sky frame (docs/sources/spec/conventions.md): origin at the telescope, z up, x north, y west.
Directions are zenith and azimuth in degrees. Zenith is measured from +z. Azimuth is
astronomical: 0 deg is +x (north) and 90 deg is -y (east), increasing clockwise seen from above.
"""

from __future__ import annotations

import math

Vector = tuple[float, float, float]


def unit_vector(zenith_deg: float, azimuth_deg: float) -> Vector:
    """Unit vector of a direction given by zenith and astronomical azimuth [deg]."""
    if not (math.isfinite(zenith_deg) and 0.0 <= zenith_deg <= 180.0):
        raise ValueError(f"zenith {zenith_deg} deg must lie in [0, 180]")
    if not math.isfinite(azimuth_deg):
        raise ValueError(f"azimuth {azimuth_deg} deg must be finite")
    zenith = math.radians(zenith_deg)
    azimuth = math.radians(azimuth_deg)

    def component(value: float) -> float:
        # Round-off below 1e-15 (sin 180 deg) is zero; adding 0.0 turns -0.0 into 0.0.
        return 0.0 if abs(value) < 1e-15 else value + 0.0

    return (
        component(math.sin(zenith) * math.cos(azimuth)),
        component(-math.sin(zenith) * math.sin(azimuth)),
        component(math.cos(zenith)),
    )


def positive(value: float, name: str) -> float:
    """``value`` as float, finite and > 0."""
    if not (math.isfinite(value) and value > 0.0):
        raise ValueError(f"{name} must be finite and > 0 (got {value})")
    return float(value)


def non_negative(value: float, name: str) -> float:
    """``value`` as float, finite and >= 0."""
    if not (math.isfinite(value) and value >= 0.0):
        raise ValueError(f"{name} must be finite and >= 0 (got {value})")
    return float(value)
