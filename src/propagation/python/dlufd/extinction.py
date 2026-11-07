"""Atmospheric extinction and early thinning, identical to sim_telarray's transmission.

Arrays are float64 for the transmission lookups and float32 for photon columns, as in the
photon stream. The lookups broadcast their array arguments.
"""

from collections.abc import Mapping, Sequence

import numpy as np

from dlufd import _core


def _vectorized(function, *arguments):
    """Broadcast float64 arguments to a common shape and apply function to flat arrays."""
    arrays = np.broadcast_arrays(*(np.asarray(a, dtype=np.float64) for a in arguments))
    shape = arrays[0].shape
    flat = [np.ascontiguousarray(a).reshape(-1) for a in arrays]
    return np.asarray(function(*flat)).reshape(shape)


class TransmissionTable(_core.TransmissionTable):
    """sim_telarray atmospheric transmission table (read_trans of sim_absorb.c).

    Args:
        path: sim_telarray transmission file (``# H2= ..., H1= ...`` header).
        site_altitude_m: observation level [m a.s.l.]; must match H2 within 25 m.
    """

    def optical_depth(self, wavelength_nm, height_m):
        """Optical depth of a vertical path from the observation level to height_m."""
        return _vectorized(super().optical_depth, wavelength_nm, height_m)

    def transmission(self, wavelength_nm, zem_m, airmass):
        """sim_telarray atmospheric_transmission: exp(-od * airmass)."""
        return _vectorized(super().transmission, wavelength_nm, zem_m, airmass)

    def transmission2(self, wavelength_nm, zem_m, airmass, od_tel, od_focus):
        """sim_telarray atmospheric_transmission2 with telescope optical depths od_tel, od_focus."""
        return _vectorized(super().transmission2, wavelength_nm, zem_m, airmass, od_tel, od_focus)

    def telescope_depths(self, telescope_altitude_m, pathlen_m=0.0):
        """Telescope optical depths (od_tel, od_focus), indexed by integer wavelength 0..999 nm.

        Args:
            telescope_altitude_m: telescope altitude a.s.l. [m].
            pathlen_m: optical path inside the telescope [m] (sim_telarray pathlen_tel).
        """
        return super().telescope_depths(float(telescope_altitude_m), float(pathlen_m))


def apply_extinction(
    columns: Mapping[str, np.ndarray],
    table: TransmissionTable,
    telescope_altitude_m: float,
    envelope_wavelength: Sequence[float] | np.ndarray,
    envelope_value: Sequence[float] | np.ndarray,
    seed: int,
    event: int,
    telescope: int,
    first_index: int,
    weighted: bool,
    wavelength_range: tuple[float, float],
    pathlen_m: float = 0.0,
) -> dict[str, np.ndarray]:
    """Extinction, wavelength assignment and envelope thinning of photons in frame ground.

    Columns (float32, equal length): x, y, cx, cy, t, zem, w; optional z, cz, wavelength
    (0 = assign from 1/lambda^2 on wavelength_range; nm, integer as in sim_telarray).

    The acceptance is p = p_atm(lambda) * E(lambda), with E the envelope (linear in lambda,
    values in [0, 1]). Thinned mode (weighted=False) returns one row per surviving photon
    with w = 1. Weighted mode returns all rows with w multiplied by p. Both return the
    column p. first_index is the index of the first photon of this telescope and event.

    Args:
        columns: photon columns by name.
        table: transmission table of the site.
        telescope_altitude_m: telescope altitude a.s.l. [m].
        envelope_wavelength, envelope_value: envelope table (strictly increasing wavelength).
        seed, event, telescope, first_index: Philox key of the draws.
        weighted: weighted (True) or thinned (False) output.
        wavelength_range: (lo, hi) [nm] for wavelengths drawn from 1/lambda^2, lo >= 1.
        pathlen_m: optical path inside the telescope [m].

    Returns:
        Dict of float32 arrays, one per output column.

    Raises:
        ValueError: invalid input (e.g. |cz| inconsistent, weights above 2048, wavelength outside
            the envelope table).
    """
    float_columns = {
        name: np.ascontiguousarray(value, dtype=np.float32) for name, value in columns.items()
    }
    return _core.apply_extinction(
        float_columns,
        table,
        float(telescope_altitude_m),
        np.ascontiguousarray(envelope_wavelength, dtype=np.float64),
        np.ascontiguousarray(envelope_value, dtype=np.float64),
        int(seed),
        int(event),
        int(telescope),
        int(first_index),
        bool(weighted),
        (float(wavelength_range[0]), float(wavelength_range[1])),
        float(pathlen_m),
    )
