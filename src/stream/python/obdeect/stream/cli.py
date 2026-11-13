"""The ``obdeect`` command: ``obdeect <component> <command> [options]``.

Each installed component registers one group in the entry-point group ``obdeect.commands``
(for example ``sources = "obdeect.sources.cli:group"``). A group is a :class:`Group` of
commands. Native commands replace the ``obdeect`` process with the program installed in the
component's ``bin`` directory, so pipes run at native speed.
"""

import argparse
import os
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from importlib import metadata, resources

ENTRY_POINTS = "obdeect.commands"
# Components and their distributions, for the hint on a missing group.
DISTRIBUTIONS = {
    "stream": "obdeect-stream",
    "propagation": "obdeect-propagation",
    "sources": "obdeect-sources",
    "ray-tracing": "obdeect-ray-tracing",
}


@dataclass(frozen=True)
class Command:
    """One command: a native program in ``<package>/bin`` or a Python function."""

    help: str
    program: str | None = None
    function: Callable[[list[str]], int] | None = None


@dataclass(frozen=True)
class Group:
    """The commands of one component."""

    package: str
    help: str
    commands: dict[str, Command] = field(default_factory=dict)


def program_path(package: str, program: str) -> str:
    """Path of a native program installed in ``<package>/bin``."""
    path = resources.files(package) / "bin" / program
    if not path.is_file():
        raise FileNotFoundError(f"{program} is not installed in {package}")
    return str(path)


def find_program(package: str, program: str) -> str | None:
    """Path of a native program in ``<package>/bin``, or None if it is not installed."""
    try:
        return program_path(package, program)
    except FileNotFoundError:
        return None


def run_command(group: Group, name: str, command: Command, argv: list[str]) -> int:
    prog = f"obdeect {name}"
    if command.program is not None:
        path = program_path(group.package, command.program)
        sys.stdout.flush()
        os.execv(path, [prog, *argv])
    return command.function(argv)


def groups() -> dict[str, Group]:
    found = {}
    for entry in metadata.entry_points(group=ENTRY_POINTS):
        found[entry.name] = entry.load()
    return dict(sorted(found.items()))


def versions() -> str:
    from obdeect.stream import format_major, format_minor

    lines = []
    for name, distribution in DISTRIBUTIONS.items():
        try:
            lines.append(f"{distribution} {metadata.version(distribution)}")
        except metadata.PackageNotFoundError:
            continue
    lines.append(f"photon-stream format {format_major}.{format_minor}")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    available = groups()
    parser = argparse.ArgumentParser(
        prog="obdeect",
        description="Run an obdeect command: obdeect <component> <command> [options].",
    )
    parser.add_argument("--version", action="store_true", help="print the installed versions")
    components = parser.add_subparsers(dest="component", metavar="component")
    for name, group in available.items():
        sub = components.add_parser(name, help=group.help, description=group.help)
        commands = sub.add_subparsers(dest="command", metavar="command", required=True)
        for command_name, command in group.commands.items():
            commands.add_parser(command_name, help=command.help, add_help=False)

    if argv and argv[0] == "--version":
        print(versions())
        return 0
    if not argv or argv[0] not in available:
        if argv and argv[0] in DISTRIBUTIONS:
            print(
                f"obdeect: {argv[0]} is not installed: pip install {DISTRIBUTIONS[argv[0]]}",
                file=sys.stderr,
            )
            return 2
        parser.parse_args(argv[:1] or ["--help"])
        return 2
    group = available[argv[0]]
    if len(argv) < 2 or argv[1] not in group.commands:
        parser.parse_args(argv[:2] if len(argv) > 1 else [argv[0], "--help"])
        return 2
    name = f"{argv[0]} {argv[1]}"
    return run_command(group, name, group.commands[argv[1]], argv[2:])


group = Group(
    package="obdeect.stream",
    help="photon streams: inspect, copy, convert from CORSIKA and to sim_telarray or GrISU",
    commands={
        "cat": Command("print or concatenate photon streams", program="obdeect-stream-cat"),
        "tee": Command("copy a stream to stdout and to files", program="obdeect-stream-tee"),
        "from-corsika": Command(
            "convert CORSIKA IACT output to a photon stream", program="obdeect-stream-from-corsika"
        ),
        "to-simtel": Command(
            "convert a photon stream to CORSIKA IACT input for sim_telarray",
            program="obdeect-stream-to-simtel",
        ),
        "to-grisu": Command(
            "convert a photon stream to GrISU photon text", program="obdeect-stream-to-grisu"
        ),
    },
)


if __name__ == "__main__":
    raise SystemExit(main())
