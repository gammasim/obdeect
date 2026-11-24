Propagation Python API: invalid arguments are rejected instead of being truncated, cast or
ignored (non-integer seeds and ids, `scattering_points` without `scattering`, `earth_radius_m`
without curvature, `wavelength_range` not a pair, an `envelope` string, float64 photon
columns). `pathlen_m` is required with extinction. A failed run leaves no output file.
