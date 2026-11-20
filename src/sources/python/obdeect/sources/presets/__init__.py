"""Presets: command-line builders of run descriptions for common light sources.

Each module has ``NAME``, ``HELP``, ``add_arguments(parser)`` and ``build(args) -> dict``; the
dictionary is a run description (see :mod:`obdeect.sources.config`).
"""

from importlib import import_module

#: Preset name -> module name.
PRESETS = {
    "flasher": "flasher",
    "led": "led",
    "point-source": "point_source",
    "beam-scan": "beam_scan",
    "laser-beam": "laser_beam",
    "lidar": "lidar",
    "night-sky": "night_sky",
    "muon": "muon",
    "shower": "shower",
}


def preset(name: str):
    """The module of a preset."""
    if name not in PRESETS:
        raise ValueError(f"unknown preset '{name}' (one of {', '.join(PRESETS)})")
    return import_module(f"obdeect.sources.presets.{PRESETS[name]}")
