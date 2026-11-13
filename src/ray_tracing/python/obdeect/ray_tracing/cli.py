"""The ``obdeect`` command: compile optical models, trace sources, plot."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from obdeect.ray_tracing import __version__
from obdeect.ray_tracing.formats import ModelError
from obdeect.ray_tracing.model import STATUS_NAMES, OpticalModel
from obdeect.ray_tracing.psf import summary
from obdeect.ray_tracing.stream import TRANSMISSION_PARTS, trace_stream
from obdeect.stream.cli import Command, Group

# Columns of `obdeect ray-tracing obscurers` (see OpticalModel.obscurer_hits).
CENSUS = ("crossed", "first", "sole")


def _compile(args: argparse.Namespace) -> None:
    model = OpticalModel.from_simulation_models(
        args.simulation_models,
        args.model,
        args.version,
        configuration_seed=args.seed,
        zenith_deg=args.zenith_deg,
    )
    model.save(args.output)
    print(
        f"{args.model} {args.version}: optical model {model.sha256[:12]} written to {args.output}"
    )
    for note in model.document["notes"]:
        print(f"note: {note}")


def _trace(args: argparse.Namespace) -> None:
    model = OpticalModel.load(args.model)
    photons = _star(model, args)
    arrivals = model.trace(
        photons, seed=args.seed, apply_responses=not args.geometry_only, threads=args.threads
    )
    if args.output:
        _write(args.output, photons, arrivals)
    print(
        json.dumps(
            {"model_sha256": model.sha256, **summary(arrivals, model.sampling_area_m2)}, indent=1
        )
    )


def _star(model: OpticalModel, args: argparse.Namespace) -> dict:
    fx, fy = (math.radians(a) for a in args.field_deg)
    towards = np.array([math.tan(fx), math.tan(fy), 1.0])
    return model.star_photons(
        args.photons,
        seed=args.seed,
        towards_source=tuple(towards / np.linalg.norm(towards)),
        distance_m=args.distance_m,
        wavelength_nm=args.wavelength_nm,
    )


def _obscurers(args: argparse.Namespace) -> None:
    model = OpticalModel.load(args.model)
    hits = model.obscurer_hits(_star(model, args), seed=args.seed, threads=args.threads)
    unobscured = hits["unobscured"]
    print(f"optical model {model.sha256[:12]}: {unobscured:.0f} photons reach the focal surface")
    print(f"without obscurers, {hits['detected']:.0f} with them; per obscurer (in percent):")
    print(f"{'name':40s} {'group':20s} {'crossed':>8s} {'first':>8s} {'sole':>8s}")
    for k in np.argsort(-hits["crossed"], kind="stable"):
        shares = (100 * hits[c][k] / unobscured if unobscured else 0.0 for c in CENSUS)
        print(
            f"{hits['name'][k]:40s} {hits['group'][k]:20s} " + " ".join(f"{s:8.3f}" for s in shares)
        )


def _write(path: Path, photons: dict, arrivals: dict) -> None:
    columns = {"photon_id": photons["photon_id"], **arrivals}
    if path.suffix == ".npz":
        np.savez_compressed(path, status_names=np.array(STATUS_NAMES), **columns)
        return
    names, data = [], []
    for name, values in columns.items():
        if values.ndim == 2:
            names += [f"{name}_{axis}" for axis in "xyz"]
            data += list(values.T)
        else:
            names.append(name)
            data.append(values)
    np.savetxt(
        path, np.column_stack(data), delimiter=",", header=",".join(names), comments="", fmt="%.17g"
    )


def _trace_stream(args: argparse.Namespace) -> None:
    if (args.pointing_zenith_deg is None) != (args.pointing_azimuth_deg is None):
        raise ModelError("give both --pointing-zenith-deg and --pointing-azimuth-deg, or neither")
    pointing = None
    if args.pointing_zenith_deg is not None:
        pointing = (args.pointing_zenith_deg, args.pointing_azimuth_deg)
    model = OpticalModel.load(args.model)
    try:
        counts = trace_stream(
            model,
            args.input,
            args.output,
            seed=args.seed,
            threads=args.threads,
            telescopes=args.telescopes,
            pointing=pointing,
            telescope_transmission=args.telescope_transmission,
        )
    except (ValueError, RuntimeError) as error:
        raise ModelError(str(error)) from error
    # Standard output may carry the stream, so the summary goes to standard error.
    print(
        f"trace-stream: {counts['events']} events, {counts['photons_in']} photons in, "
        f"{counts['photons_out']} out",
        file=sys.stderr,
    )


def _telescope_ids(text: str) -> tuple[int, ...]:
    try:
        ids = tuple(int(item) for item in text.split(",") if item.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"telescope ids must be integers: {text!r}") from error
    if not ids:
        raise argparse.ArgumentTypeError("give at least one telescope id")
    return ids


def _plot(args: argparse.Namespace) -> None:
    try:
        from obdeect.ray_tracing.plot import plot
    except ModuleNotFoundError as error:
        if error.name != "matplotlib":
            raise
        raise ModelError(
            "plotting requires Matplotlib; install obdeect-ray-tracing[plot]"
        ) from error

    arrivals = dict(np.load(args.arrivals)) if args.arrivals else None
    plot(OpticalModel.load(args.model), arrivals, args.output)
    print(f"plot written to {args.output}")


def _source_arguments(command: argparse.ArgumentParser) -> None:
    command.add_argument("--model", type=Path, required=True)
    command.add_argument(
        "--seed", type=int, required=True, help="seed of source sampling and mirror scatter"
    )
    command.add_argument("--photons", type=int, default=100000)
    command.add_argument(
        "--field-deg",
        type=float,
        nargs=2,
        default=(0.0, 0.0),
        metavar=("X", "Y"),
        help="source direction: field angles from the optical axis",
    )
    command.add_argument(
        "--distance-m", type=float, default=math.inf, help="source distance (default: star)"
    )
    command.add_argument("--wavelength-nm", type=float, default=400.0)
    command.add_argument("--threads", type=int, default=0, help="0 uses every core")


def main(argv: list[str] | None = None, prog: str = "obdeect ray-tracing") -> int:
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    c = commands.add_parser("compile", help="compile an optical model from simulation-models")
    c.add_argument(
        "--simulation-models", type=Path, required=True, help="simulation-models checkout"
    )
    c.add_argument("--model", required=True, help="telescope model, e.g. LSTN-design")
    c.add_argument("--version", required=True, help="production version, e.g. 7.0.0")
    c.add_argument("--seed", type=int, required=True, help="seed of the random mirror errors")
    c.add_argument(
        "--zenith-deg", type=float, help="telescope zenith angle for zenith-dependent errors"
    )
    c.add_argument("--output", type=Path, required=True)
    c.set_defaults(run=_compile)

    t = commands.add_parser("trace", help="trace a star or point source")
    _source_arguments(t)
    t.add_argument("--geometry-only", action="store_true", help="do not apply optical responses")
    t.add_argument("--output", type=Path, help="per-photon results (.npz or .csv)")
    t.set_defaults(run=_trace)

    o = commands.add_parser(
        "obscurers",
        help="which obscurers shade a star: the share of the light each one meets (crossed), "
        "absorbs (first) or alone stops (sole)",
    )
    _source_arguments(o)
    o.set_defaults(run=_obscurers)

    s = commands.add_parser(
        "trace-stream",
        help="trace a photon stream (G or T photons) to camera-frame (C) photons",
    )
    s.add_argument("--model", type=Path, required=True)
    s.add_argument("--seed", type=int, required=True, help="run seed of the random draws")
    s.add_argument("--threads", type=int, default=0, help="0 uses every core")
    s.add_argument(
        "--telescopes",
        type=_telescope_ids,
        metavar="IDS",
        help="comma-separated telescope ids to trace; other telescopes are dropped",
    )
    s.add_argument(
        "--pointing-zenith-deg",
        type=float,
        help="telescope zenith [deg]; default: the run keys telescope_pointing_zenith_rad "
        "and telescope_pointing_azimuth_rad",
    )
    s.add_argument(
        "--pointing-azimuth-deg",
        type=float,
        help="telescope azimuth [deg] in the array frame (x north, y west), as simtools' "
        "telescope phi; not the CORSIKA shower azimuth; default: the run keys",
    )
    s.add_argument(
        "--telescope-transmission",
        choices=TRANSMISSION_PARTS,
        default="full",
        help="full: apply the telescope transmission, including its on-axis value; angular: "
        "only t(theta)/t0, leaving t0 pending in the ledger for the camera",
    )
    s.add_argument("-i", "--input", default="-", help="input stream (default: standard input)")
    s.add_argument("-o", "--output", default="-", help="output stream (default: standard output)")
    s.set_defaults(run=_trace_stream)

    p = commands.add_parser("plot", help="plot the dish layout and the focal-plane image")
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--arrivals", type=Path, help="results of 'obdeect ray-tracing trace' (.npz)")
    p.add_argument("--output", type=Path, required=True)
    p.set_defaults(run=_plot)

    args = parser.parse_args(argv)
    try:
        args.run(args)
    except (ModelError, OSError) as error:
        print(f"obdeect: error: {error}", file=sys.stderr)
        return 1
    return 0


def _command(name: str):
    return lambda argv: main([name, *argv])


group = Group(
    package="obdeect.ray_tracing",
    help="ray tracing through telescope optics",
    commands={
        "compile": Command(
            "compile an optical model from simulation-models", function=_command("compile")
        ),
        "trace": Command("trace a star or point source", function=_command("trace")),
        "obscurers": Command("which obscurers shade a star", function=_command("obscurers")),
        "trace-stream": Command(
            "trace a photon stream to camera-frame photons", function=_command("trace-stream")
        ),
        "plot": Command(
            "plot the dish layout and the focal-plane image", function=_command("plot")
        ),
    },
)
