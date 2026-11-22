"""Atmosphere stage on a photon stream: extinction, envelope thinning, refraction.

The stream is read from and written to files; the work is done by the C++ core, the same code
as the ``obdeect propagation run`` command-line tool.
"""

import operator
import os
from collections.abc import Iterable, Mapping, Sequence

from obdeect.propagation import _core

MODES = ("thinned", "weighted")
CURVATURES = ("none", "corsika", "spherical")
WAVELENGTH_INTERPOLATIONS = ("linear", "floor")
DEFAULT_SCATTERING_POINTS = 8


def propagate(
    input_path: str | os.PathLike,
    output_path: str | os.PathLike,
    *,
    transmission: str | os.PathLike | None = None,
    site_altitude_m: float,
    seed: int,
    profile: str | os.PathLike | None = None,
    envelope: Mapping[str, str | os.PathLike] | Sequence[tuple[str, str | os.PathLike]] = (),
    wavelength_range: tuple[float, float] | None = None,
    mode: str | None = None,
    store_p: bool = False,
    scattering: str | os.PathLike | None = None,
    scattering_points: int | None = None,
    scattering_primaries: Iterable[int] | None = None,
    pathlen_m: float | None = None,
    telescope_pathlen_m: Mapping[int, float] | None = None,
    extinction: bool = True,
    curvature: str = "none",
    earth_radius_m: float | None = None,
    telescopes: Iterable[int] | None = None,
    threads: int = 1,
    wavelength_interpolation: str = "linear",
    dispersion_reference_nm: float | None = None,
) -> dict[str, float]:
    """Run the atmosphere stage on a photon stream file.

    Args:
        input_path: input photon stream.
        output_path: output photon stream (overwritten).
        transmission: transmission table of the site (extinction-table format, ``# H2=`` header).
            Required with extinction; not allowed without.
        site_altitude_m: observation level [m a.s.l.]; must match the table within 25 m.
        seed: Philox seed of the wavelength and thinning draws.
        profile: optional atmospheric profile (atmprof format); enables the air delay and
            refraction of sky-frame photons (unless the run says refraction is applied).
        envelope: envelope factors, as a mapping name -> file or (name, file) pairs. Each file
            is a two-column text table (wavelength nm, value). The product is applied.
        wavelength_range: (lo, hi) nm for wavelengths drawn from 1/lambda^2. Default: the run
            keys wavelength_min_nm and wavelength_max_nm.
        mode: ``"thinned"`` (one row per surviving photon, w = 1) or ``"weighted"``
            (w multiplied by the acceptance). Default ``"thinned"``; without extinction the
            output is always weighted and ``"thinned"`` is an error.
        store_p: keep the acceptance column p in the output.
        scattering: optional scattering table (docs/scattering.md). Enables single scattering of
            the sky-frame photons before extinction; scattered photons go to their telescope's
            block. The run needs telescope_r_m. Off when None.
        scattering_points: scattering points per photon (1 to 1024); used only with scattering.
            Default None: 8 points.
        scattering_primaries: telescope ids (1-based) whose photons are scattered; needs
            scattering. Default None: the photons of all telescopes. Ids must exist in the run.
        pathlen_m: optical path inside each telescope [m] (the in-telescope part of the
            extinction, see docs/propagation.md). Default None: 0 m.
        telescope_pathlen_m: per-telescope path [m] as {telescope_id: metres}, ids 1-based;
            overrides pathlen_m for those telescopes. Default None.
        extinction: False for geometry only: sky to ground frame (refraction, air delay) and
            wavelength assignment, no extinction, envelope or thinning. All photons keep their
            weight and the efficiency status records ``atmosphere`` as pending. Cannot be
            combined with transmission, envelope, mode="thinned", scattering, store_p or the
            path options.
        curvature: Earth curvature of the atmosphere. ``"none"`` (default):
            plane-parallel, input with ``curvature_applied = 1`` is rejected. ``"corsika"``:
            photons from CORSIKA with the CURVED option (``curvature_applied = 1``, frame G);
            only the curved slant depth replaces the plane-parallel airmass (CORSIKA applied the
            rest). ``"spherical"``: photons
            without curvature (``curvature_applied = 0``); the curvature is applied to the
            emission altitudes and the local vertical of each telescope; with a profile, the
            refraction and air delay follow the spherical geometry. Not combined with
            scattering.
        earth_radius_m: Earth radius [m] for the curved modes; must be above 1e6 m. Default
            None: the CORSIKA CURVED radius (6371315 m).
        telescopes: telescope ids (1-based) to keep. The blocks (telescope, photons,
            telescope_end and their efficiency status) of other telescopes are dropped from the
            output, and their photons are not decoded. The run-level records are unchanged.
            With scattering, the kept telescopes are the scattering targets; other blocks are
            still read for primaries but not written. Counts in the result are for the kept
            telescopes. Default None: all telescopes.
        threads: worker threads, 1 to 256. The output is identical for any count. Default 1.
        wavelength_interpolation: how the transmission table is read at a wavelength:
            ``"linear"`` (default, continuous in lambda: the depths and the envelope are
            interpolated between the integer nm rows) or ``"floor"`` (the row at floor(lambda), as
            sim_telarray). Needs transmission.
        dispersion_reference_nm: wavelength [nm] (200 to 2000) of the profile's n-1. Scales the
            bending and the travel time of the bent sky-frame photons to their wavelength
            (docs/propagation/refraction.md, "Dispersion"). Needs profile; rejected for a run
            with refraction_applied = 1. Default None: no dispersion.

    Returns:
        Dict with photons_in, photons_out, dropped (sky-frame photons not reaching the
        telescope plane), dropped_by_reason (dict of the dropped photons by reason: moving_away,
        horizontal, below_ground, turns_over, not_bendable, misses_level; the values sum to
        dropped) and weight_out (sum of w of the output). See docs/propagation/propagation.md.

    Raises:
        ValueError: invalid settings or input, or an input efficiency status that shows a
            factor as applied.
        RuntimeError: unreadable or truncated files.
    """
    if mode is not None and mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, not {mode!r}")
    if curvature not in CURVATURES:
        raise ValueError(f"curvature must be one of {CURVATURES}, not {curvature!r}")
    if scattering_primaries is not None and scattering is None:
        raise ValueError("scattering_primaries needs scattering")
    if extinction:
        if transmission is None:
            raise ValueError("transmission is required (or extinction=False)")
    else:
        _check_geometry_only(
            transmission,
            envelope,
            mode,
            scattering,
            scattering_points,
            store_p,
            pathlen_m,
            telescope_pathlen_m,
        )
    items = envelope.items() if isinstance(envelope, Mapping) else envelope
    factors = [(str(name), os.fspath(path)) for name, path in items]
    # operator.index rejects non-integer ids (1.5, "1") instead of truncating them.
    pathlens = {
        operator.index(telescope): float(metres)
        for telescope, metres in (telescope_pathlen_m or {}).items()
    }
    selection = None
    if telescopes is not None:
        selection = [operator.index(telescope) for telescope in telescopes]
        if not selection:
            raise ValueError("telescopes: give at least one telescope id")
    threads = operator.index(threads)
    if not 1 <= threads <= 256:
        raise ValueError("threads must be an integer from 1 to 256")
    if wavelength_interpolation not in WAVELENGTH_INTERPOLATIONS:
        raise ValueError(
            f"wavelength_interpolation must be one of {WAVELENGTH_INTERPOLATIONS}, "
            f"not {wavelength_interpolation!r}"
        )
    if dispersion_reference_nm is not None and profile is None:
        raise ValueError("dispersion_reference_nm needs profile")
    if not extinction and wavelength_interpolation != "linear":
        raise ValueError("extinction=False cannot be combined with wavelength_interpolation")
    return dict(
        _core.propagate(
            os.fspath(input_path),
            os.fspath(output_path),
            None if transmission is None else os.fspath(transmission),
            float(site_altitude_m),
            None if profile is None else os.fspath(profile),
            factors,
            None
            if wavelength_range is None
            else (float(wavelength_range[0]), float(wavelength_range[1])),
            int(seed),
            mode == "weighted" or not extinction,
            bool(store_p),
            None if scattering is None else os.fspath(scattering),
            DEFAULT_SCATTERING_POINTS if scattering_points is None else int(scattering_points),
            0.0 if pathlen_m is None else float(pathlen_m),
            pathlens,
            not extinction,
            curvature,
            None if earth_radius_m is None else float(earth_radius_m),
            selection,
            threads,
            None
            if scattering_primaries is None
            else [operator.index(telescope) for telescope in scattering_primaries],
            wavelength_interpolation,
            None if dispersion_reference_nm is None else float(dispersion_reference_nm),
        )
    )


def to_sky(
    input_path: str | os.PathLike,
    output_path: str | os.PathLike,
    *,
    profile: str | os.PathLike | None = None,
) -> dict[str, int]:
    """Convert the frame G photons of a photon stream to frame sky (the inverse of the stage).

    The command ``obdeect propagation to-sky`` with the same arguments. For a run with a profile
    (refraction_applied = 1), ``propagate(..., extinction=False, profile=...)`` on the output
    reproduces the input's x, y, cx, cy, t and zem to float32 rounding; see
    docs/propagation/ground-to-sky.md.

    Args:
        input_path: input photon stream with frame G photons (for example from CORSIKA).
        output_path: output photon stream in frame sky (overwritten).
        profile: the run's atmospheric profile; inverts the air delay. Needs
            refraction_applied = 1 in the run (photons already refracted); None: no air delay.

    Returns:
        Dict with photons_in, photons_out and dropped (rows that cannot be converted).

    Raises:
        ValueError: invalid input (not frame G, curved run, bending that is not inverted).
        RuntimeError: unreadable or truncated files.
    """
    from obdeect.propagation import _core

    return dict(
        _core.to_sky(
            os.fspath(input_path),
            os.fspath(output_path),
            None if profile is None else os.fspath(profile),
        )
    )


def _check_geometry_only(
    transmission,
    envelope,
    mode,
    scattering,
    scattering_points,
    store_p,
    pathlen_m,
    telescope_pathlen_m,
) -> None:
    """Reject settings that need extinction or thinning (extinction=False)."""
    rejected = {
        "transmission": transmission is not None,
        "envelope": bool(envelope),
        'mode="thinned"': mode == "thinned",
        "scattering": scattering is not None,
        "scattering_points": scattering_points is not None,
        "store_p": bool(store_p),
        "pathlen_m": pathlen_m is not None,
        "telescope_pathlen_m": bool(telescope_pathlen_m),
    }
    given = [name for name, is_set in rejected.items() if is_set]
    if given:
        raise ValueError(f"extinction=False cannot be combined with {', '.join(given)}")
