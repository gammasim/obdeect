"""laser-beam: a pulsed laser pointing up or down at a position above the observation level.

Model: one point emission of the laser pulse, ``E lambda / (h c)`` primary photons with pulse
energy ``E`` (``--pulse-energy``) and wavelength ``--wavelength-nm``. It starts at
``--laser-position`` (or at ``--start-height`` on the beam), with the beam axis, the divergence
(full angle) and the pulse profile of the options. Scattering is computed by propagation.

Chain::

    obdeect sources run laser-beam ... -o primaries.pst
    obdeect propagation run --transmission T --site-altitude H --seed N --mode weighted \\
        --scattering S --scattering-points 64 -i primaries.pst -o scattered.pst

t = 0 is the laser firing. Limits: a non-horizontal beam's tilt from the vertical (nadir for a
downward beam, the zenith for an upward one) is at most ``--max-tilt``; a horizontal beam has no
tilt limit and starts ``--start-distance`` along the beam (see ``_pulsed_beam``). The atmosphere is
plane-parallel and only single scattering is computed, so the scattering chain is for beams within
``--max-tilt``.
Direction convention and the direct-beam check: see ``_pulsed_beam``.
"""

from __future__ import annotations

import argparse

from obdeect.sources.presets import _common, _pulsed_beam

NAME = "laser-beam"
HELP = "pulsed laser (up or down): primary photons of one pulse (point emission)"


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the laser-beam preset (no device defaults: all geometry is required)."""
    _common.add_run_arguments(parser)
    group = parser.add_argument_group("laser")
    group.add_argument(
        "--laser-position",
        type=_common.vector,
        required=True,
        metavar="X,Y,Z",
        help="laser position in the sky frame; z above the observation level [m]",
    )
    _pulsed_beam.add_arguments(parser)


def build(args: argparse.Namespace) -> dict:
    """Run description of one pulse per event."""
    return _pulsed_beam.build(args, position=args.laser_position)
