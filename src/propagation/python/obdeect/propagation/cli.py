"""Commands of ``obdeect propagation``."""

from obdeect.stream.cli import Command, Group


def _convert_table(argv: list[str]) -> int:
    from obdeect.propagation.convert import main

    return main(argv)


group = Group(
    package="obdeect.propagation",
    help="atmospheric propagation: extinction, refraction, timing, thinning and scattering",
    commands={
        "run": Command(
            "propagate a photon stream through the atmosphere", program="obdeect-propagation-run"
        ),
        "convert-table": Command(
            "convert a transmission or scattering table", function=_convert_table
        ),
    },
)
