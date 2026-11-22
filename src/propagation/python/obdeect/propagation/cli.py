"""Commands of ``obdeect propagation``."""

from obdeect.stream.cli import Command, Group


def _convert_table(argv: list[str]) -> int:
    from obdeect.propagation.convert import main

    return main(argv)


def _pathlen(argv: list[str]) -> int:
    from obdeect.propagation.optics_path import main

    return main(argv)


group = Group(
    package="obdeect.propagation",
    help="atmospheric propagation: extinction, refraction, timing, thinning and scattering",
    commands={
        "run": Command(
            "propagate a photon stream through the atmosphere", program="obdeect-propagation-run"
        ),
        "to-sky": Command(
            "convert ground-frame photons to the sky frame", program="obdeect-propagation-to-sky"
        ),
        "convert-table": Command(
            "convert a transmission or scattering table", function=_convert_table
        ),
        "pathlen": Command(
            "print the in-telescope path length [m] of an optical model", function=_pathlen
        ),
    },
)
