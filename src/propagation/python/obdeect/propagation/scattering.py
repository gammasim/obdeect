"""Single scattering of light-source photons (Rayleigh and Henyey-Greenstein).

See docs/propagation/scattering.md.

Arrays are float64 for the table lookups and float32 for photon columns, as in the photon stream.
"""

from collections.abc import Mapping, Sequence

import numpy as np
from obdeect.propagation import _core
from obdeect.propagation._arrays import broadcast_call, float32_columns, integer


class ScatteringTable(_core.ScatteringTable):
    """Scattering table: molecular and aerosol optical depths, aerosol albedo and asymmetry.

    Args:
        path: scattering table file (columns: wavelength_nm altitude_km tau_mol tau_aer ssa g).
    """

    def profile(self, wavelength_nm, altitude_m) -> tuple[np.ndarray, ...]:
        """Return (tau_mol, tau_aer, ssa, g) at wavelength [nm] and altitude [m a.s.l.]."""
        return broadcast_call(super().profile, wavelength_nm, altitude_m)


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

    Each primary photon is scattered at n_points points along its path (sampled by optical depth):
    downward to the observation level, or upward (cz > 0) to the top of the scattering table;
    horizontal photons (cz = 0) are not scattered. For each point and telescope one scattered
    photon is aimed at the sphere (through a uniform point of its projected disc). It travels down
    if the point is above the sphere centre and up if it is below. Its weight is
    w * P_scatter * phase * Omega_tel * T(point -> telescope). Position and zem are those of the
    scattering point; t is t_in + (emission -> point) / c in vacuum, without the path from the point
    to the telescope plane (the atmosphere stage adds it).

    Columns (float32, equal length): x, y, z (z above the observation level, m), cx, cy, t, w,
    wavelength (integer nm, > 0); optional cz (either sign) and zem (consistency checked). Other
    dtypes are rejected. The primary photons are not changed; their extinction is applied by the
    extinction step. Scattered photons carry T in their weight: do not apply the atmospheric
    acceptance to them.

    Args:
        columns: primary photon columns by name.
        table: scattering table.
        transmission: obdeect.propagation.TransmissionTable of the site (gives the observation level
            and T).
        telescopes: (telescope_id, x_m, y_m, z_m, radius_m) per telescope sphere; z above the
            observation level, z_m >= 0.
        seed, event: Philox key of the draws (integers).
        first_index: index of the first primary photon of the event (integer).
        n_points: scattering points per primary photon (1 to 1024, integer).

    Returns:
        Dict keyed by telescope_id of dicts of float32 arrays (x, y, z, cx, cy, cz, t, zem, w,
        wavelength), one row per scattered photon.
    """
    sphere = [
        (integer(i, "telescope_id"), float(x), float(y), float(z), float(r))
        for i, x, y, z, r in telescopes
    ]
    return _core.scatter(
        float32_columns(columns),
        table,
        transmission,
        sphere,
        integer(seed, "seed"),
        integer(event, "event"),
        integer(first_index, "first_index"),
        integer(n_points, "n_points"),
    )
