"""Single scattering of light-source photons (Rayleigh and Henyey-Greenstein).

See docs/scattering.md.

Arrays are float64 for the table lookups and float32 for photon columns, as in the photon stream.
"""

from collections.abc import Mapping, Sequence

import numpy as np

from dlufd import _core


def _vectorized(function, *arguments):
    """Broadcast float64 arguments to a common shape and apply function to flat arrays."""
    arrays = np.broadcast_arrays(*(np.asarray(a, dtype=np.float64) for a in arguments))
    shape = arrays[0].shape
    flat = [np.ascontiguousarray(a).reshape(-1) for a in arrays]
    return tuple(np.asarray(r).reshape(shape) for r in function(*flat))


class ScatteringTable(_core.ScatteringTable):
    """Scattering table: molecular and aerosol optical depths, aerosol albedo and asymmetry.

    Args:
        path: scattering table file (columns: wavelength_nm altitude_km tau_mol tau_aer ssa g).
    """

    def profile(self, wavelength_nm, altitude_m):
        """Return (tau_mol, tau_aer, ssa, g) at wavelength [nm] and altitude [m a.s.l.]."""
        return _vectorized(super().profile, wavelength_nm, altitude_m)


def scatter(
    columns: Mapping[str, np.ndarray],
    table: ScatteringTable,
    transmission,
    telescopes: Sequence[tuple[int, float, float, float, float]],
    seed: int,
    event: int,
    first_index: int,
    n_points: int,
) -> dict[int, dict[str, np.ndarray]]:
    """Single scattering of primary photons (frame sky) towards telescope spheres.

    Each primary photon is scattered at n_points points along its downward path (sampled by
    optical depth). For each point and telescope one scattered photon is aimed at the sphere
    centre (through a uniform point of its projected disc). Its weight is
    w * P_scatter * phase * Omega_tel * T(point -> telescope). Position and zem are those of the
    scattering point; t is t_in + (emission -> point) / c in vacuum, without the path from the point
    to the telescope plane (the atmosphere stage adds it).

    Columns (float32, equal length): x, y, z (z above the observation level, m), cx, cy, t, w,
    wavelength (integer nm, > 0); optional cz and zem (consistency checked).
    The primary photons are not changed; their extinction is applied by the extinction step.
    Scattered photons carry T in their weight: do not apply the atmospheric acceptance to them.

    Args:
        columns: primary photon columns by name.
        table: scattering table.
        transmission: dlufd.TransmissionTable of the site (gives the observation level and T).
        telescopes: (telescope_id, x_m, y_m, z_m, radius_m) per telescope sphere; z above the
            observation level, z_m >= 0.
        seed, event: Philox key of the draws.
        first_index: index of the first primary photon of the event.
        n_points: scattering points per primary photon (1 to 1024).

    Returns:
        Dict keyed by telescope_id of dicts of float32 arrays (x, y, z, cx, cy, cz, t, zem, w,
        wavelength), one row per scattered photon.
    """
    float_columns = {
        name: np.ascontiguousarray(value, dtype=np.float32) for name, value in columns.items()
    }
    sphere = [(int(i), float(x), float(y), float(z), float(r)) for i, x, y, z, r in telescopes]
    return _core.scatter(
        float_columns,
        table,
        transmission,
        sphere,
        int(seed),
        int(event),
        int(first_index),
        int(n_points),
    )
