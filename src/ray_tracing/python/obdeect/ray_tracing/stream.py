"""Trace a photon stream through an optical model (the Python API of ``trace-stream``).

The tracing is the native code of the program ``obdeect-ray-tracing-trace-stream``: this module
hands the optical models and the file paths to ``_core.trace_stream``, which releases the GIL while
it runs. The command and the function therefore write the same bytes. Ground-frame (G) or
telescope-frame (T) photons go in; camera-frame (C) photons come out. Records pass through
unchanged, except that telescopes not selected are dropped and the efficiency status is updated
(docs/ray-tracing/photon-stream.md).
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
from obdeect.ray_tracing import _core
from obdeect.ray_tracing.errors import ModelError
from obdeect.ray_tracing.model import OpticalModel

# Efficiency status factor left for a downstream module when only the angular part is applied.
CONSTANT_FACTOR = "telescope_transmission_constant"


def _model_text(model: OpticalModel) -> str:
    """The JSON text of an optical model, as OpticalModel.save writes it."""
    return json.dumps(model.document, ensure_ascii=False, indent=1) + "\n"


def _model_sources(
    model: OpticalModel | Mapping[int, OpticalModel],
) -> list[tuple[int | None, str]]:
    """(telescope id or None for every telescope, model text) of each model of a run."""
    if isinstance(model, OpticalModel):
        return [(None, _model_text(model))]
    items = [(int(telescope), value) for telescope, value in model.items()]
    if not items:
        raise ModelError("give at least one optical model")
    # A model shared by several telescopes is serialised once.
    texts = {id(value): _model_text(value) for _, value in items}
    return [(telescope, texts[id(value)]) for telescope, value in items]


def _stream_path(path: Any, standard: str) -> str:
    """The file of a stream; "-" is the process's standard input or output."""
    text = str(path)
    return standard if text == "-" else text


def _input_planes(model: OpticalModel, theta_deg: float, phi_deg: float) -> tuple[float, float]:
    try:
        return _core.input_planes(_model_text(model), float(theta_deg), float(phi_deg))
    except ValueError as error:
        raise ModelError(str(error)) from error


def ground_to_telescope_rotation(theta_deg: float, phi_deg: float) -> np.ndarray:
    """Rotation R with p_telescope = R p_ground for a telescope at zenith theta, azimuth phi.

    The ground frame has x north, y west and z up (the array frame, as simtools
    ``ground_direction_to_telescope_frame``). phi is the telescope azimuth in that
    frame, not the CORSIKA propagation azimuth of a shower. The telescope frame has
    +z along the optical axis. A shower that propagates along -axis maps to
    direction (0, 0, -1).
    """
    return np.asarray(_core.ground_to_telescope_rotation(float(theta_deg), float(phi_deg)))


def minimum_input_plane_height_m(
    model: OpticalModel, theta_deg: float, phi_deg: float = 0.0
) -> float:
    """Height [m] of the ground-frame input plane above the telescope sphere centre.

    The plane lies above every incoming obscurer and every mirror in the pointing direction,
    with the axes offset and a margin included. A photon starts on the plane and moves down to
    the ground, so a mirror point above the plane would never be reached: at grazing pointing
    the dish reaches far above the obscurers.
    """
    return _input_planes(model, theta_deg, phi_deg)[0]


def lower_input_plane_height_m(
    model: OpticalModel, theta_deg: float, phi_deg: float = 0.0
) -> float:
    """Height [m] of the ground-frame input plane for photons from below.

    The plane lies a margin below the lowest optic (any obscurer, the primary and the
    secondary) along the zenith of the pointing, with the axes offset included. Photons
    from below start on it and move up; it may lie below the ground.
    """
    return _input_planes(model, theta_deg, phi_deg)[1]


def trace_stream(
    model: OpticalModel | Mapping[int, OpticalModel],
    input: Any,
    output: Any,
    *,
    seed: int,
    threads: int = 0,
    telescopes: Iterable[int] | None = None,
    pointing: tuple[float, float] | None = None,
    telescope_transmission: str = "full",
) -> dict[str, int]:
    """Trace the photons of a photon stream and write the camera-frame photons.

    Parameters
    ----------
    model: OpticalModel or mapping of int to OpticalModel
        Optical model of every traced telescope, or the optical model of each telescope id.
        With a mapping, telescopes without a model are dropped, and the efficiency status of
        each telescope block names the model that traced it.
    input, output: str or path-like
        Photon-stream files; "-" is standard input or output.
    seed: int
        Run seed (non-negative). Draws are keyed on (seed, event, telescope, photon), as in the
        other components.
    threads: int
        Worker threads (0 uses every core); results do not depend on it.
    telescopes: iterable of int or None
        Telescope ids to trace; other telescopes are dropped. None traces all (with a mapping
        of models: all telescopes with a model). Each id needs a model.
    pointing: (float, float) or None
        Telescope (zenith, azimuth) [deg]; azimuth in the array frame (x north,
        y west), as simtools' telescope phi. None takes it from the run keys
        telescope_pointing_zenith_rad and telescope_pointing_azimuth_rad [rad].
        Ground-frame photons need one of the two; the event header has no pointing.
    telescope_transmission: "full" or "angular"
        "full" applies the telescope transmission t(theta), including its on-axis
        value t0. "angular" applies t(theta) / t0 only; the efficiency status then leaves
        t0 as the pending factor ``telescope_transmission_constant`` for the
        downstream camera module. Models without telescope transmission reject
        "angular".

    Returns
    -------
    dict
        Counts of ``events``, ``photons_in`` and ``photons_out``.

    Raises
    ------
    ModelError
        For an invalid stream, a missing pointing (ground-frame photons only), an
        unsupported frame or column, an efficiency status that says a factor is already
        applied, an unknown telescope_transmission, "angular" with a model that has none, and
        a telescope id without a model or outside 1..65535.
    """
    if seed < 0:
        raise ModelError("seed must not be negative")
    selected = None if telescopes is None else [int(telescope) for telescope in telescopes]
    point = None if pointing is None else (float(pointing[0]), float(pointing[1]))
    try:
        events, photons_in, photons_out = _core.trace_stream(
            _stream_path(input, "/dev/stdin"),
            _stream_path(output, "/dev/stdout"),
            _model_sources(model),
            seed=int(seed),
            threads=int(threads),
            telescopes=selected,
            pointing=point,
            telescope_transmission=telescope_transmission,
        )
    except (ValueError, RuntimeError) as error:
        raise ModelError(str(error)) from error
    return {"events": events, "photons_in": photons_in, "photons_out": photons_out}
