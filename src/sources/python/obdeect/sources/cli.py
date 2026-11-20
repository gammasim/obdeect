"""Commands of ``obdeect sources``.

``run SOURCE`` writes the photon stream of a run; ``check SOURCE`` validates a run and prints its
expected counts. SOURCE is a ``.toml`` or ``.json`` run file, or the name of a preset (its options
follow SOURCE). ``--set KEY=VALUE`` overrides one value of the run description, for example
``--set events[0].repeat=5``.
"""

import argparse
import copy
import importlib
import json
import re
import sys
from collections.abc import Sequence
from typing import Any

from obdeect.sources import config, presets
from obdeect.stream.cli import Command, Group

RUN_SUFFIXES = (".toml", ".json")
_TOKEN = re.compile(r"\.?([^.\[\]]+)|\[(\d+)\]")


class _UsageError(Exception):
    """A command-line error, reported on one line."""


class _Parser(argparse.ArgumentParser):
    def error(self, message: str):
        raise _UsageError(message)


def parse_set(text: str) -> tuple[list[str | int], Any]:
    """``KEY=VALUE`` of ``--set``: a path such as ``events[0].repeat`` and the value.

    The value is parsed as JSON (numbers, true, lists, ...); anything that is not valid JSON is a
    string.
    """
    key, eq, raw = text.partition("=")
    if not eq or not key:
        raise ValueError(f"--set '{text}': expected KEY=VALUE")
    try:
        value = json.loads(raw)
    except ValueError:
        value = raw
    return _path(key, text), value


def _path(key: str, text: str) -> list[str | int]:
    tokens: list[str | int] = []
    pos = 0
    while pos < len(key):
        match = _TOKEN.match(key, pos)
        if match is None or match.end() == pos:
            raise ValueError(f"--set '{text}': bad path at '{key[pos:]}'")
        tokens.append(match.group(1) if match.group(1) is not None else int(match.group(2)))
        pos = match.end()
    if not tokens:
        raise ValueError(f"--set '{text}': empty path")
    return tokens


def _item(node: Any, token: str | int, text: str) -> Any:
    """The child ``token`` of ``node``; a missing key of a table is created as an empty table."""
    if isinstance(token, int):
        if not isinstance(node, list):
            raise ValueError(f"--set '{text}': [{token}] needs a list")
        if token >= len(node):
            raise ValueError(f"--set '{text}': index {token} out of range (length {len(node)})")
        return node[token]
    if not isinstance(node, dict):
        raise ValueError(f"--set '{text}': '{token}' needs a table")
    return node.setdefault(token, {})


def apply_set(description: dict, path: list[str | int], value: Any, text: str) -> None:
    """Set the value at ``path`` in ``description`` (changed in place)."""
    node = description
    for token in path[:-1]:
        node = _item(node, token, text)
    last = path[-1]
    if isinstance(last, int):
        _item(node, last, text)  # checks the list and the index
        node[last] = value
    else:
        if not isinstance(node, dict):
            raise ValueError(f"--set '{text}': '{last}' needs a table")
        node[last] = value


def _overrides(description: dict, args: argparse.Namespace) -> dict:
    """Copy of the run description with ``--seed`` and then ``--set`` applied."""
    result = copy.deepcopy(description)
    seed = getattr(args, "seed_override", None)
    if seed is not None:
        apply_set(result, ["run", "seed"], seed, f"--seed {seed}")
    for text in getattr(args, "overrides", None) or []:
        path, value = parse_set(text)
        apply_set(result, path, value, text)
    return result


def _is_file(source: str) -> bool:
    return source.lower().endswith(RUN_SUFFIXES)


def _preset_names() -> str:
    return ", ".join(presets.PRESETS)


def _add_common(parser: argparse.ArgumentParser, command: str) -> None:
    """Options of ``run`` and ``check``; an option that a preset defines stays the preset's."""
    options = [
        (
            ("--set",),
            {
                "dest": "overrides",
                "action": "append",
                "default": [],
                "metavar": "KEY=VALUE",
                "help": "override a value of the run description, e.g. "
                "events[0].repeat=5 (value as JSON, else a string); repeatable",
            },
        )
    ]
    if command == "run":
        options.append((
            ("-o", "--output"),
            {
                "metavar": "OUT.pst",
                "default": None,
                "help": "photon-stream output file (default and '-': stdout)",
            },
        ))
        options.append((
            ("--max-rows",),
            {
                "type": float,
                "default": 1e8,
                "metavar": "N",
                "help": "refuse runs expected to write more photon rows (all telescopes; "
                "40 bytes per row; default 1e8)",
            },
        ))
        options.append((
            ("--seed",),
            {
                "dest": "seed_override",
                "type": int,
                "metavar": "N",
                "help": "random seed of the run (shortcut for --set run.seed=N)",
            },
        ))
    else:
        options.append((
            ("--json",),
            {
                "action": "store_true",
                "help": "print the normalized run description as JSON (usable as a run file)",
            },
        ))
    for flags, kwargs in options:
        try:
            parser.add_argument(*flags, **kwargs)
        except argparse.ArgumentError:
            pass  # a preset defines this option


def _top_parser(command: str) -> argparse.ArgumentParser:
    parser = _Parser(
        prog=f"obdeect sources {command}",
        usage=f"obdeect sources {command} SOURCE [options]",
        description=f"{command} SOURCE: a run file (.toml or .json), or a preset.",
        epilog="presets:\n"
        + "\n".join(f"  {name:<14}{presets.preset(name).HELP}" for name in presets.PRESETS)
        + f"\n\nRun 'obdeect sources {command} PRESET --help' for the options of a preset.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common(parser, command)
    return parser


def _source_parser(command: str, source: str, module: Any | None) -> argparse.ArgumentParser:
    parser = _Parser(
        prog=f"obdeect sources {command} {source}",
        description=getattr(module, "HELP", None),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    if module is not None:
        module.add_arguments(parser)
    _add_common(parser, command)
    return parser


def _load(command: str, argv: list[str]) -> tuple[dict, argparse.Namespace]:
    """Run description of SOURCE (file or preset) with the overrides, and the parsed options."""
    source, rest = argv[0], argv[1:]
    if source.startswith("-"):
        raise _UsageError(f"SOURCE must come first, got '{source}'")
    if _is_file(source):
        parser = _source_parser(command, source, None)
        args = parser.parse_args(rest)
        description = config.load(source)
    else:
        module = presets.preset(source)
        parser = _source_parser(command, source, module)
        args = parser.parse_args(rest)
        description = module.build(args)
    return _overrides(description, args), args


def _print_check(result: dict) -> None:
    """Counts and expected weights: scalar values first, then one row per emission."""
    emissions = result.get("emissions", [])
    for key, value in result.items():
        if key != "emissions":
            print(f"{key}: {value}")
    if not emissions:
        return
    columns = list(emissions[0])
    rows = [[str(item.get(column, "")) for column in columns] for item in emissions]
    widths = [max(len(c), *(len(r[i]) for r in rows)) for i, c in enumerate(columns)]
    print("  ".join(c.ljust(w) for c, w in zip(columns, widths)))
    for row in rows:
        print("  ".join(v.ljust(w) for v, w in zip(row, widths)))


def _execute(command: str, argv: list[str]) -> int:
    if not argv:
        raise _UsageError(
            f"missing SOURCE (a .toml or .json run file, or a preset: {_preset_names()})"
        )
    if argv[0] in ("-h", "--help"):
        _top_parser(command).print_help()
        return 0
    description, args = _load(command, argv)
    if command == "check":
        if args.json:
            print(json.dumps(config.normalize(description), indent=2))
            return 0
        module = importlib.import_module("obdeect.sources.run")
        _print_check(module.check(description))
        return 0

    module = importlib.import_module("obdeect.sources.run")
    output = "-" if args.output is None else args.output
    sys.stdout.flush()
    summary = module.write(description, output, max_rows=args.max_rows)
    sys.stdout.flush()
    print(
        f"events {summary['events']}, rows {summary['rows']}, "
        f"sum of weights {summary['weight']:.6g}, expected {summary['expected']:.6g}",
        file=sys.stderr,
    )
    return 0


def _command(name: str, argv: Sequence[str]) -> int:
    try:
        return _execute(name, list(argv))
    except (ValueError, OSError, _UsageError) as error:
        print(f"obdeect sources {name}: {error}", file=sys.stderr)
        return 1


def run(argv: list[str]) -> int:
    """``obdeect sources run``."""
    return _command("run", argv)


def check(argv: list[str]) -> int:
    """``obdeect sources check``."""
    return _command("check", argv)


group = Group(
    package="obdeect.sources",
    help="light sources: run descriptions from presets or run files, and their photon streams",
    commands={
        "run": Command("write the photon stream of a run (run file or preset)", function=run),
        "check": Command("validate a run and print its expected counts", function=check),
    },
)
