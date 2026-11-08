"""Atmosphere stage on a gstrickd photon stream: extinction, envelope thinning, refraction.

The stream is read from and written to files; the work is done by the C++ core, the same code
as the ``dlufd-propagate`` command-line tool.
"""

import operator
import os
from collections.abc import Mapping, Sequence

from dlufd import _core

MODES = ("thinned", "weighted")


def propagate(
    input_path: str | os.PathLike,
    output_path: str | os.PathLike,
    *,
    transmission: str | os.PathLike,
    site_altitude_m: float,
    seed: int,
    profile: str | os.PathLike | None = None,
    envelope: Mapping[str, str | os.PathLike] | Sequence[tuple[str, str | os.PathLike]] = (),
    wavelength_range: tuple[float, float] | None = None,
    mode: str = "thinned",
    store_p: bool = False,
    scattering: str | os.PathLike | None = None,
    scattering_points: int = 8,
    pathlen_m: float = 0.0,
    telescope_pathlen_m: Mapping[int, float] | None = None,
) -> dict[str, float]:
    """Run the atmosphere stage on a photon stream file.

    Args:
        input_path: input gstrickd stream.
        output_path: output gstrickd stream (overwritten).
        transmission: transmission table of the site (extinction-table format, ``# H2=`` header).
        site_altitude_m: observation level [m a.s.l.]; must match the table within 25 m.
        seed: Philox seed of the wavelength and thinning draws.
        profile: optional atmospheric profile (atmprof format); enables the air delay and
            refraction of sky-frame photons (unless the run says refraction is applied).
        envelope: envelope factors, as a mapping name -> file or (name, file) pairs. Each file
            is a two-column text table (wavelength nm, value). The product is applied.
        wavelength_range: (lo, hi) nm for wavelengths drawn from 1/lambda^2. Default: the run
            keys wavelength_min_nm and wavelength_max_nm.
        mode: ``"thinned"`` (one row per surviving photon, w = 1) or ``"weighted"``
            (w multiplied by the acceptance).
        store_p: keep the acceptance column p in the output.
        scattering: optional scattering table (docs/scattering.md). Enables single scattering of
            the sky photons before extinction; scattered photons go to their telescope's block.
            The run needs telescope_r_m. Off when None.
        scattering_points: scattering points per photon (1 to 1024); used only with scattering.
        pathlen_m: optical path inside each telescope [m] (the in-telescope part of the
            extinction, see docs/propagation.md). Default 0.
        telescope_pathlen_m: per-telescope path [m] as {telescope_id: metres}, ids 1-based;
            overrides pathlen_m for those telescopes. Default None.

    Returns:
        Dict with photons_in, photons_out, dropped (sky photons not reaching the telescope
        plane) and weight_out (sum of w of the output).

    Raises:
        ValueError: invalid settings or input, or an input ledger that shows a factor as applied.
        RuntimeError: unreadable or truncated files.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, not {mode!r}")
    items = envelope.items() if isinstance(envelope, Mapping) else envelope
    factors = [(str(name), os.fspath(path)) for name, path in items]
    # operator.index rejects non-integer ids (1.5, "1") instead of truncating them.
    pathlens = {
        operator.index(telescope): float(metres)
        for telescope, metres in (telescope_pathlen_m or {}).items()
    }
    return dict(
        _core.propagate(
            os.fspath(input_path),
            os.fspath(output_path),
            os.fspath(transmission),
            float(site_altitude_m),
            None if profile is None else os.fspath(profile),
            factors,
            None
            if wavelength_range is None
            else (float(wavelength_range[0]), float(wavelength_range[1])),
            int(seed),
            mode == "weighted",
            bool(store_p),
            None if scattering is None else os.fspath(scattering),
            int(scattering_points),
            float(pathlen_m),
            pathlens,
        )
    )
