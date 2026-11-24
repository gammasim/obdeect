"""Atmosphere stage on a photon stream: extinction, envelope thinning, refraction.

The stream is read from and written to files; the work is done by the C++ core, the same code
as the ``obdeect propagation run`` command-line tool. The core validates the settings; this module
checks what only the Python call can show (a setting that is given but would be ignored) and
converts the arguments.
"""

import os
from collections.abc import Iterable, Mapping, Sequence

from obdeect.propagation import _core
from obdeect.propagation._arrays import float_pair, integer

#: Default scattering points per photon (used only with scattering).
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
    ground_representation: str = "plane",
) -> dict[str, float]:
    """Run the atmosphere stage on a photon stream file.

    Args:
        input_path: input photon stream.
        output_path: output photon stream (overwritten).
        transmission: transmission table of the site (extinction-table format, ``# H2=`` header).
            Required with extinction; not allowed without.
        site_altitude_m: observation level [m a.s.l.]; must match the table within 25 m.
        seed: Philox seed of the wavelength and thinning draws (integer).
        profile: optional atmospheric profile (atmprof format); enables the air delay and
            refraction of sky-frame photons (unless the run says refraction is applied).
        envelope: envelope factors, as a mapping name -> file or a sequence of (name, file) pairs.
            Each file is a two-column text table (wavelength nm, value). The product is applied.
        wavelength_range: (lo, hi) nm for wavelengths drawn from 1/lambda^2. Default: the run
            keys wavelength_min_nm and wavelength_max_nm.
        mode: ``"thinned"`` (one row per surviving photon, w = 1) or ``"weighted"``
            (w multiplied by the acceptance). Default ``"thinned"``; without extinction the
            output is always weighted and ``"thinned"`` is an error.
        store_p: keep the acceptance column p in the output.
        scattering: optional scattering table (docs/propagation/scattering.md). Enables single
            scattering of the sky-frame photons before extinction; scattered photons go to their
            telescope's block. The run needs telescope_r_m. Off when None.
        scattering_points: scattering points per photon (1 to 1024); needs scattering. Default
            None: 8 points.
        scattering_primaries: telescope ids (1-based) whose photons are scattered; needs
            scattering. Default None: the photons of all telescopes. Ids must exist in the run.
        pathlen_m: optical path inside each telescope [m] (the in-telescope part of the
            extinction, see docs/propagation/extinction.md). Required with extinction (0 for
            none); rejected without extinction. See ``obdeect.propagation.in_telescope_path_m``.
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
            rest). ``"spherical"``: photons without curvature (``curvature_applied = 0``); the
            curvature is applied to the emission altitudes and the local vertical of each
            telescope; with a profile, the refraction and air delay follow the spherical geometry.
            Not combined with scattering.
        earth_radius_m: Earth radius [m] for the curved modes; must be above 1e6 m. Default
            None: the CORSIKA CURVED radius (6371315 m). Rejected with curvature ``"none"``.
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
        ground_representation: ``"plane"`` (default: rows at z = 0 relative to the telescope) or
            ``"sphere"`` (rows where the path enters the telescope sphere, columns z, cz, zem and
            path; docs/propagation/propagation.md, "Sphere representation"). Sphere needs the run
            key telescope_r_m. Sphere rows take extinction along the chord through the sphere
            (docs/propagation/extinction.md, "Chord integral").

    Returns:
        Dict with photons_in, photons_out, dropped (sky-frame photons not reaching the
        telescope plane), dropped_by_reason (dict of the dropped photons by reason: moving_away,
        horizontal, below_ground, turns_over, not_bendable, misses_level, misses_sphere; the
        values sum to dropped), thinned_out (weight removed by the acceptance draw of the
        primaries), scattered_dropped and scattered_dropped_by_reason (the scattered photons that do
        not reach the telescope, by reason), scattered_thinned_out (weight of the scattered photons
        removed by the acceptance) and weight_out (sum of w of the output). See
        docs/propagation/propagation.md, "Accounting".

    Raises:
        ValueError: invalid settings or input, a setting that is missing or given where it is not
            used, or an input efficiency status that shows a factor as applied.
        TypeError: an id or the seed is not an integer, or a path is not a path.
        RuntimeError: unreadable or truncated files.
    """
    if extinction:
        if transmission is None:
            raise ValueError("transmission is required (or extinction=False)")
        if pathlen_m is None:
            raise ValueError(
                "pathlen_m is required with extinction: the in-telescope path [m] "
                "(0 for none; see obdeect.propagation.in_telescope_path_m)"
            )
    else:
        # Settings the core cannot see: their defaults would hide that they were given.
        rejected = {
            "scattering_points": scattering_points is not None,
            "pathlen_m": pathlen_m is not None,
            "telescope_pathlen_m": bool(telescope_pathlen_m),
            "wavelength_interpolation": wavelength_interpolation != "linear",
        }
        given = [name for name, is_set in rejected.items() if is_set]
        if given:
            raise ValueError(f"extinction=False cannot be combined with {', '.join(given)}")
    if earth_radius_m is not None and curvature == "none":
        raise ValueError("earth_radius_m needs curvature corsika or spherical")
    if scattering_points is not None and scattering is None:
        raise ValueError("scattering_points needs scattering")

    # The mode and the other settings are validated by the core (the same messages as the
    # command-line tool); a geometry-only run defaults to weighted output.
    core_mode = mode if mode is not None else ("thinned" if extinction else "weighted")
    pathlens = {
        integer(telescope, "telescope_pathlen_m key"): float(metres)
        for telescope, metres in (telescope_pathlen_m or {}).items()
    }
    return dict(
        _core.propagate(
            os.fspath(input_path),
            os.fspath(output_path),
            None if transmission is None else os.fspath(transmission),
            float(site_altitude_m),
            None if profile is None else os.fspath(profile),
            _envelope_pairs(envelope),
            None if wavelength_range is None else float_pair(wavelength_range, "wavelength_range"),
            integer(seed, "seed"),
            core_mode,
            bool(store_p),
            None if scattering is None else os.fspath(scattering),
            DEFAULT_SCATTERING_POINTS
            if scattering_points is None
            else integer(scattering_points, "scattering_points"),
            0.0 if pathlen_m is None else float(pathlen_m),
            pathlens,
            not extinction,
            curvature,
            None if earth_radius_m is None else float(earth_radius_m),
            None if telescopes is None else [integer(t, "telescopes") for t in telescopes],
            integer(threads, "threads"),
            None
            if scattering_primaries is None
            else [integer(t, "scattering_primaries") for t in scattering_primaries],
            wavelength_interpolation,
            None if dispersion_reference_nm is None else float(dispersion_reference_nm),
            ground_representation,
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
    return dict(
        _core.to_sky(
            os.fspath(input_path),
            os.fspath(output_path),
            None if profile is None else os.fspath(profile),
        )
    )


def _envelope_pairs(envelope) -> list[tuple[str, str]]:
    """The envelope factors as (name, file) pairs of strings. A mapping or a sequence of pairs."""
    if isinstance(envelope, str | bytes):
        raise ValueError("envelope must be a mapping or a sequence of (name, file) pairs")
    items = list(envelope.items()) if isinstance(envelope, Mapping) else list(envelope)
    pairs = []
    for item in items:
        if isinstance(item, str | bytes):
            raise ValueError(f"envelope entry must be a (name, file) pair, not {item!r}")
        try:
            name, path = item
        except TypeError, ValueError:
            raise ValueError(f"envelope entry must be a (name, file) pair, not {item!r}") from None
        if not isinstance(name, str):
            raise TypeError(f"envelope factor name must be a string, not {name!r}")
        pairs.append((name, os.fspath(path)))
    return pairs
