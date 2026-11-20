"""Beam and field-of-view geometry of the beam presets (laser-beam, lidar).

Sky frame (docs/sources/spec/conventions.md): origin at the telescope, z up, x north, y west.
Directions are zenith and azimuth in degrees. Zenith is measured from +z. Azimuth is
astronomical: 0 deg is +x (north) and 90 deg is -y (east), increasing clockwise seen from above.

A telescope field of view is a cone with its apex at the origin, its axis along the pointing
direction and a half opening angle. Only the forward nappe (axial coordinate s > 0) counts.
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
    # Adding 0.0 turns a negative zero into 0.0 (cosmetic, for the run description).
    return (
        math.sin(zenith) * math.cos(azimuth) + 0.0,
        -math.sin(zenith) * math.sin(azimuth) + 0.0,
        math.cos(zenith) + 0.0,
    )


def half_angle(fov_deg: float) -> float:
    """Half opening angle [deg] of a full field-of-view angle in (0, 90) deg."""
    if not (math.isfinite(fov_deg) and 0.0 < fov_deg < 90.0):
        raise ValueError(f"field of view {fov_deg} deg must lie in (0, 90)")
    return 0.5 * fov_deg


def _dot(u, v) -> float:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _roots(a: float, b: float, c: float) -> tuple[float, float] | None:
    """Sorted real roots of a t^2 + b t + c with a != 0, or None without real roots.

    Uses q = -(b + sign(b) sqrt(disc)) / 2 so that neither root suffers from cancellation.
    """
    disc = b * b - 4.0 * a * c
    if disc < 0.0:
        return None
    q = -0.5 * (b + math.copysign(math.sqrt(disc), b))
    if q == 0.0:
        return (0.0, 0.0)
    r1, r2 = q / a, c / q
    return (min(r1, r2), max(r1, r2))


def cone_path(
    origin: Vector,
    direction: Vector,
    axis: Vector,
    half_angle_deg: float,
    t_max: float = math.inf,
) -> tuple[float, float]:
    """Entry distance and path length of the ray ``origin + t * direction`` inside a cone.

    The cone apex is the origin of the frame, ``axis`` and ``direction`` are unit vectors, and
    t runs over [0, t_max]. Returns ``(t_enter, length)``.

    Inside the cone means: axial coordinate ``s > 0`` and ``|perpendicular| <= tan(half) s``.
    Writing the ray as ``p + t d`` gives ``a t^2 + b t + c <= 0`` with
    ``a = |dp|^2 - T^2 dz^2``, ``b = 2 (pp.dp - T^2 pz dz)``, ``c = |pp|^2 - T^2 pz^2``
    (``T = tan(half)``, subscript z axial, p perpendicular). Its roots are computed without
    cancellation. The cone is convex, so the forward intersection is a single interval.

    Raises ValueError if the ray misses the cone (no interval of positive length) or if the
    intersection is unbounded (the ray stays inside the cone for all t up to infinity).
    """
    if not (math.isfinite(half_angle_deg) and 0.0 < half_angle_deg < 90.0):
        raise ValueError("half opening angle must lie in (0, 90) deg")
    miss = "beam does not cross the field of view"
    unbounded = "beam has an unbounded intersection with the field of view"

    tan2 = math.tan(math.radians(half_angle_deg)) ** 2
    dz = _dot(direction, axis)
    pz = _dot(origin, axis)
    dp = tuple(d - dz * x for d, x in zip(direction, axis))
    pp = tuple(p - pz * x for p, x in zip(origin, axis))
    a = _dot(dp, dp) - tan2 * dz * dz
    b = 2.0 * (_dot(pp, dp) - tan2 * pz * dz)
    c = _dot(pp, pp) - tan2 * pz * pz

    # Forward nappe: pz + t dz > 0, and 0 <= t <= t_max.
    lo, hi = 0.0, t_max
    if dz > 0.0:
        lo = max(lo, -pz / dz)
    elif dz < 0.0:
        hi = min(hi, -pz / dz)
    elif pz <= 0.0:
        raise ValueError(miss)

    # a vanishes up to rounding when the ray is parallel to a cone generator.
    a_scale = max(1.0, _dot(dp, dp), tan2 * dz * dz)
    if abs(a) > 1e-12 * a_scale:
        roots = _roots(a, b, c)
        if a > 0.0:
            if roots is None:
                raise ValueError(miss)
            lo, hi = max(lo, roots[0]), min(hi, roots[1])
        elif roots is not None:
            # Inside only outside the roots; the forward part is the single non-empty piece.
            first = (lo, min(hi, roots[0]))
            second = (max(lo, roots[1]), hi)
            lo, hi = first if first[1] > first[0] else second
    elif b > 0.0:
        hi = min(hi, -c / b)
    elif b < 0.0:
        lo = max(lo, -c / b)
    elif c > 0.0:
        raise ValueError(miss)

    if hi <= lo:
        raise ValueError(miss)
    if math.isinf(hi):
        raise ValueError(unbounded)
    return lo, hi - lo


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
