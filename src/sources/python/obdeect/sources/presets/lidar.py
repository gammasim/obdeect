"""lidar: one satellite lidar pulse seen by a ground telescope (one line emission per pulse).

The satellite is at ``--satellite-position`` (or at ``(0, 0, --satellite-altitude)``). Its z is
the height above the observation level, and the telescope is at the origin. The beam leaves the
satellite along ``--beam-zenith``/``--beam-azimuth``, which must point towards the ground
(zenith > 90 deg). Its full length to the ground is ``L = z_sat / (-cz)``.

The emission is one line from the satellite. Its path is clipped to the telescope
field-of-view cone and to the ground (0 <= t <= L). Its photons are the photons per pulse
times ``clipped length / L``. Photons per pulse are ``E lambda / (h c)`` with the pulse energy
``E`` and the wavelength ``lambda``. The beam has full divergence ``--divergence`` (filled cone
with half angle divergence/2; cone 0 for zero). Only primary photons are generated; scattering
happens in obdeect propagation. ``--events`` is the number of pulses (one event per pulse).

Azimuth convention: astronomical, clockwise from north seen from above (sky frame x north, y
west), as in laser-beam.
"""

from __future__ import annotations

import argparse

from obdeect.sources.presets import _common
from obdeect.sources.presets._beam import (
    cone_path,
    half_angle,
    non_negative,
    positive,
    unit_vector,
)

NAME = "lidar"
HELP = "satellite lidar pulse seen by a ground telescope (one line emission per pulse)"

#: Planck constant [J s] and speed of light [m/s] (exact SI values).
PLANCK_J_S = 6.62607015e-34
SPEED_OF_LIGHT_M_S = 299792458.0


def photons_per_pulse(energy_j: float, wavelength_nm: float) -> float:
    """Photons in a pulse of energy ``E`` [J] at wavelength ``lambda`` [nm]: E lambda / (h c)."""
    return energy_j * (wavelength_nm * 1e-9) / (PLANCK_J_S * SPEED_OF_LIGHT_M_S)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the lidar preset (no device defaults: all geometry is required)."""
    _common.add_run_arguments(parser)
    group = parser.add_argument_group("satellite")
    where = group.add_mutually_exclusive_group(required=True)
    where.add_argument(
        "--satellite-position",
        type=_common.vector,
        metavar="X,Y,Z",
        help="satellite position in the sky frame; z above the observation level [m]",
    )
    where.add_argument(
        "--satellite-altitude",
        type=float,
        metavar="M",
        help="satellite height above the observation level, directly above the telescope [m]",
    )
    group.add_argument(
        "--beam-zenith",
        type=float,
        required=True,
        metavar="DEG",
        help="zenith of the beam propagation direction; > 90 (towards the ground) [deg]",
    )
    group.add_argument(
        "--beam-azimuth",
        type=float,
        required=True,
        metavar="DEG",
        help="astronomical azimuth of the beam: clockwise from north [deg]",
    )
    group.add_argument(
        "--divergence",
        type=float,
        required=True,
        metavar="DEG",
        help="full beam divergence angle, in [0, 180) [deg]",
    )
    group = parser.add_argument_group("pulse")
    group.add_argument(
        "--pulse-energy", type=float, required=True, metavar="J", help="pulse energy [J]"
    )
    group.add_argument(
        "--wavelength-nm",
        type=float,
        required=True,
        metavar="NM",
        help="wavelength [nm]; the spectrum is fixed at this wavelength",
    )
    group.add_argument(
        "--pulse",
        default="fixed:time_ns=0",
        metavar="KIND[:K=V,...]",
        help="time profile of the pulse (default fixed:time_ns=0)",
    )
    group.add_argument(
        "--rows",
        type=float,
        default=1e6,
        metavar="N",
        help="sampled photon rows per pulse (a sampling choice; weights carry the physical "
        "photons; default 1e6)",
    )
    group = parser.add_argument_group("telescope")
    group.add_argument(
        "--telescope-zenith",
        type=float,
        required=True,
        metavar="DEG",
        help="zenith of the telescope pointing (cone axis) [deg]",
    )
    group.add_argument(
        "--telescope-azimuth",
        type=float,
        required=True,
        metavar="DEG",
        help="astronomical azimuth of the telescope pointing: clockwise from north [deg]",
    )
    group.add_argument(
        "--telescope-radius", type=float, required=True, metavar="M", help="telescope radius [m]"
    )
    group.add_argument(
        "--fov-deg",
        type=float,
        required=True,
        metavar="DEG",
        help="full opening angle of the field of view, in (0, 90) [deg]",
    )


def build(args: argparse.Namespace) -> dict:
    """Run description with one line emission per event (one event per pulse)."""
    half = half_angle(args.fov_deg)
    radius = positive(args.telescope_radius, "telescope radius")
    energy = positive(args.pulse_energy, "pulse energy")
    wavelength = positive(args.wavelength_nm, "wavelength")
    rows = positive(args.rows, "rows")
    divergence = non_negative(args.divergence, "divergence")
    if divergence >= 180.0:
        raise ValueError(f"divergence {divergence} deg must be < 180")

    if args.satellite_position is not None:
        satellite = tuple(args.satellite_position)
    else:
        satellite = (0.0, 0.0, positive(args.satellite_altitude, "satellite altitude"))
    if satellite[2] <= 0.0:
        raise ValueError(f"satellite z {satellite[2]} m must be above the observation level")

    beam = unit_vector(args.beam_zenith, args.beam_azimuth)
    if beam[2] >= 0.0:
        raise ValueError("beam must point towards the ground: --beam-zenith must be > 90 deg")
    axis = unit_vector(args.telescope_zenith, args.telescope_azimuth)
    ground = satellite[2] / -beam[2]  # full length from the satellite to z = 0
    enter, length = cone_path(satellite, beam, axis, half, t_max=ground)
    start = [satellite[i] + beam[i] * enter for i in range(3)]

    per_pulse = photons_per_pulse(energy, wavelength)
    photons = per_pulse * length / ground
    if divergence == 0.0:
        angular = {"kind": "cone", "angle_deg": 0.0}
    else:
        angular = {"kind": "filled_cone", "angle_deg": 0.5 * divergence}
    emission = {
        "kind": "line",
        "start_m": start,
        "direction": list(beam),
        "length_m": length,
        "beta": 1.0,
        "angular": angular,
        "photons": photons,
        "bunch_size": max(1.0, photons / rows),
        "spectrum": {"kind": "fixed", "wavelength_nm": wavelength},
        "pulse": _common.profile(args.pulse, "pulse"),
        "target": "none",
    }
    return {
        "run": _common.run_section(args, frame="sky"),
        "telescopes": [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": radius}],
        "events": [{"repeat": args.events, "emissions": [emission]}],
    }
