"""Commands of ``obdeect sources``."""

from obdeect.stream.cli import Command, Group

PROGRAMS = {
    "flasher": "camera flasher (dish or secondary mirror)",
    "led": "LEDs inside the camera",
    "point-source": "point source at a finite distance",
    "beam-scan": "beam scanned across the sky",
    "laser-beam": "laser beam seen from the side",
    "night-sky": "night-sky background and stars",
    "muon-track": "Cherenkov light of a muon track",
}

group = Group(
    package="obdeect.sources",
    help="light sources: calibration devices, stars, night sky and muons",
    commands={
        name: Command(help, program=f"obdeect-sources-{name}") for name, help in PROGRAMS.items()
    },
)
