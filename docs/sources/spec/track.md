# Track source

Header `obdeect/sources/track.hpp`. A line segment emitting Cherenkov light at a fixed angle to its
direction of motion; only photons that reach a sphere are generated (importance sampling).

`Track{start (upper end), direction (unit), length_m, cos_theta (= 1/(n beta)), beta, time_ns
(emission time at start), spectrum, photons (mean emitted into 4 pi), target SphereTarget}`.

## Validation (`validate(track)`)

`photons` finite, `0..2^53`; start within limits, direction finite and unit to 1e-9;
`0 <= length` within limits; `0 < cos_theta <= 1`; `0 <= beta <= 1`; `time_ns` within limits;
sphere centre within limits and radius positive-within-limits; spectrum is Fixed
(positive-within-limits) or InverseSquare (`0 < low < high`, within limits), nothing else.

## Hit arc

For emission point `P = start + direction s`: `D = |centre - P|`; if `D <= R` every azimuth hits
(`hit`, `centre 0`, `half pi`). Else `cos_alpha = sqrt(1 - (R/D)^2)`, `unit = (centre - P)/D`,
`cos_b = dot(direction, unit)`, `perp = unit - direction cos_b`, `sin_b = |perp|`, `phi0 =
atan2(dot(perp, e2), dot(perp, e1))` with `(e1,e2) = orthonormal_basis(direction)`,
`a = cos_alpha - cos_theta cos_b`, `b = sin_theta sin_b`. If `b <= 0`: all azimuths if `a < 0`,
none otherwise. Else `x = a/b`: none if `x >= 1`; all if `x <= -1`; else `centre = phi0`,
`half = acos(x)`. Photons with azimuth within `phi0 +- half` hit; fraction `half/pi`.
Test: against brute-force azimuth scan within 2e-5.

## Emission (`emit_track(track, seed, event, telescope, out)`)

Validates; `n = track_count` (rounding rule with mean `photons`). For photon `k = 0..n-1` (all
draws keyed by k): `s = length * target_u0`; compute the arc; if no hit, skip (no row, index
`k` is consumed). Else `phi = centre + half (2 direction_u0 - 1)`; position `start + direction
s`; direction `direction cos_theta + (e1 cos phi + e2 sin phi) sin_theta`; time `time_ns + s /
(beta c)` if `beta > 0` else `time_ns`; weight `half/pi`; wavelength from the spectrum (wavelength
stream). The emission angle to the axis is `acos(cos_theta)` to 1e-12.

`expected_track_weight(track, steps = 20000)` = `photons * mean over midpoints of half/pi`
(0 for zero length or zero photons) (midpoint rule over `s`).
