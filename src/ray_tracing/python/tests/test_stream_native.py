"""The command ``obdeect ray-tracing trace-stream`` and ``trace_stream`` write the same bytes.

Both run the same native code (obdeect-ray-tracing-trace-stream; the Python function calls it
through _core). Each case traces one input both ways and compares the output files byte for byte,
or checks that both reject the input with the same message. The command-line options are tested
here too.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from obdeect.ray_tracing import OpticalModel
from obdeect.ray_tracing.formats import ModelError
from obdeect.ray_tracing.stream import trace_stream
from telescopes import dual, segmented
from test_stream import (
    RUN,
    SEED,
    beam,
    parallel_beam_in_both_frames,
    with_direction,
    write_stream,
)


def obdeect_command() -> str | None:
    """The obdeect command of this interpreter's environment (its trace-stream runs the native
    program), or None if it is not installed there."""
    beside = Path(sys.executable).with_name("obdeect")
    if beside.exists():
        return str(beside)
    return shutil.which("obdeect")


@pytest.fixture(scope="module")
def command():
    path = obdeect_command()
    if path is None:
        pytest.skip("the obdeect command is not installed")
    return path


@pytest.fixture(scope="module")
def models(tmp_path_factory):
    folder = tmp_path_factory.mktemp("models")
    out = {}
    for name, build in (("segmented", segmented), ("dual", dual)):
        records, files = build()
        model = OpticalModel.from_records(
            records, files.__getitem__, configuration_seed=1, propagation_group_index=1.0
        )
        path = folder / f"{name}.json"
        model.save(path)
        out[name] = (OpticalModel.load(path), path)
    return out


def command_argv(
    command,
    model_args,
    source,
    output,
    *,
    telescopes=None,
    pointing=(0.0, 0.0),
    transmission="full",
    threads=1,
):
    """Arguments of `obdeect ray-tracing trace-stream` for one trace."""
    argv = [command, "ray-tracing", "trace-stream", "--seed", str(SEED), "--threads", str(threads)]
    for item in model_args:
        argv += ["--model", item]
    if telescopes is not None:
        argv += ["--telescopes", ",".join(map(str, telescopes))]
    if pointing is not None:
        argv += ["--pointing-zenith-deg", repr(pointing[0])]
        argv += ["--pointing-azimuth-deg", repr(pointing[1])]
    return argv + ["--telescope-transmission", transmission, "-i", str(source), "-o", str(output)]


def run_command(command, model_args, source, output, **options):
    """Run the command; returns the completed process (not checked)."""
    argv = command_argv(command, model_args, source, output, **options)
    return subprocess.run(argv, capture_output=True, text=True, check=False)


def compare(
    command,
    models,
    tmp_path,
    source,
    *,
    model="segmented",
    telescopes=None,
    pointing=(0.0, 0.0),
    transmission="full",
    threads=1,
):
    """Trace `source` both ways: equal output bytes, or the same error from both.

    Returns the counts of the Python function, or the error message.
    """
    if isinstance(model, dict):
        python_model = {tel: models[name][0] for tel, name in model.items()}
        model_args = [f"{tel}={models[name][1]}" for tel, name in model.items()]
    else:
        python_model = models[model][0]
        model_args = [str(models[model][1])]
    python_out, command_out = tmp_path / "python.pst", tmp_path / "command.pst"
    try:
        python_counts = trace_stream(
            python_model,
            source,
            python_out,
            seed=SEED,
            threads=threads,
            telescopes=telescopes,
            pointing=pointing,
            telescope_transmission=transmission,
        )
        python_error = None
    except (ModelError, ValueError, RuntimeError) as error:  # RuntimeError: malformed stream input
        python_error = str(error)
    result = run_command(
        command,
        model_args,
        source,
        command_out,
        telescopes=telescopes,
        pointing=pointing,
        transmission=transmission,
        threads=threads,
    )
    if python_error is not None:
        assert result.returncode == 1, result.stderr
        assert python_error in result.stderr
        return python_error
    assert result.returncode == 0, result.stderr
    summary = (
        f"trace-stream: {python_counts['events']} events, {python_counts['photons_in']} photons "
        f"in, {python_counts['photons_out']} out"
    )
    assert summary in result.stderr
    assert python_out.read_bytes() == command_out.read_bytes()
    return python_counts


def test_ground_frame_on_and_off_axis(command, models, tmp_path):
    photons = beam(4000, (0.02, -0.01), radius_m=6.0)
    photons["zem"] = np.linspace(8000.0, 12000.0, 4000, dtype=np.float32)
    source = write_stream(tmp_path / "in.pst", {1: photons}, blob=True)
    counts = compare(command, models, tmp_path, source)
    assert counts["photons_out"] > 0


def test_mixed_directions_match(command, models, tmp_path):
    down = with_direction(beam(3000, (0.0, 0.0), radius_m=0.5), 0.0, 0.0, -1.0)
    up = with_direction(beam(1500, (0.01, 0.0), radius_m=0.5, seed=8), 0.01, 0.0, +1.0)
    photons = {name: np.concatenate([down[name], up[name]]) for name in down}
    source = write_stream(tmp_path / "in.pst", {1: photons})
    counts = compare(command, models, tmp_path, source)
    assert counts["photons_in"] == 4500
    assert counts["photons_out"] > 0


def test_horizontal_photons_error_matches(command, models, tmp_path):
    photons = with_direction(beam(5, (1.0, 0.0)), 1.0, 0.0, 1.0)
    photons["cz"] = np.zeros(5, np.float32)
    source = write_stream(tmp_path / "in.pst", {1: photons})
    error = compare(command, models, tmp_path, source)
    assert isinstance(error, str) and "horizontal photons" in error


def test_stdin_matches_file_input(command, models, tmp_path):
    source = write_stream(tmp_path / "in.pst", {1: beam(4000, (0.0, 0.0), radius_m=6.0)})
    model_args = [str(models["segmented"][1])]
    options = dict(telescopes=None, pointing=(0.0, 0.0), transmission="full", threads=1)
    file_result = subprocess.run(
        command_argv(command, model_args, source, "-", **options),
        capture_output=True,
        check=False,
    )
    stdin_result = subprocess.run(
        command_argv(command, model_args, "-", "-", **options),
        input=source.read_bytes(),
        capture_output=True,
        check=False,
    )
    assert file_result.returncode == 0, file_result.stderr
    assert stdin_result.returncode == 0, stdin_result.stderr
    assert b"4000 photons in" in stdin_result.stderr
    assert stdin_result.stderr == file_result.stderr
    assert stdin_result.stdout == file_result.stdout


@pytest.mark.parametrize("pointing", [(12.0, 73.0), (35.0, 180.0)])
def test_tilted_pointing(command, models, tmp_path, pointing):
    for frame, photons in zip(
        (2, 1), parallel_beam_in_both_frames(models["segmented"][0], *pointing, 3000, 5.0)
    ):
        source = write_stream(tmp_path / f"in{frame}.pst", {1: photons}, frame=frame)
        compare(command, models, tmp_path, source, pointing=pointing)


def test_run_pointing_and_dual_mirror(command, models, tmp_path):
    run = dict(RUN, telescope_pointing_zenith_rad=0.3, telescope_pointing_azimuth_rad=1.1)
    photons = parallel_beam_in_both_frames(
        models["dual"][0], 0.3 * 180 / np.pi, 1.1 * 180 / np.pi, 3000, 2.0
    )[1]
    source = write_stream(tmp_path / "in.pst", {1: photons}, run=run)
    compare(command, models, tmp_path, source, model="dual", pointing=None)


def test_envelope_status_and_other_factors(command, models, tmp_path):
    qe = {
        "factors": "qe,mirror_reflectivity",
        "qe.state": "envelope",
        "qe.wavelength_nm": np.array([300.0, 600.0]),
        "qe.value": np.array([0.3, 0.4]),
        "qe.source": "sim_telarray",
        "mirror_reflectivity.state": "envelope",
        "mirror_reflectivity.wavelength_nm": np.array([300.0, 450.0, 600.0]),
        "mirror_reflectivity.value": np.array([1.5, 2.0, 3.0]),
    }
    photons = beam(5000, (0.0, 0.0), radius_m=6.0)
    photons["wavelength"] = np.linspace(310.0, 590.0, 5000, dtype=np.float32)
    tel = {"factors": "telescope_transmission", "telescope_transmission.state": "pending"}
    source = write_stream(
        tmp_path / "in.pst",
        {1: photons},
        run_efficiency_statuses=[qe],
        telescope_efficiency_statuses={1: tel},
    )
    compare(command, models, tmp_path, source)


@pytest.mark.parametrize("transmission", ["full", "angular"])
def test_transmission_parts(command, models, tmp_path, transmission):
    source = write_stream(tmp_path / "in.pst", {1: beam(3000, (0.01, 0.0), radius_m=6.0)})
    compare(command, models, tmp_path, source, transmission=transmission)


def test_per_telescope_models_and_selection(command, models, tmp_path):
    photons = {
        1: beam(2000, (0.0, 0.0), radius_m=6.0),
        2: beam(2000, (0.0, 0.0), radius_m=2.0),
        3: beam(50, (0.0, 0.0)),
    }
    source = write_stream(tmp_path / "in.pst", photons)
    compare(
        command,
        models,
        tmp_path,
        source,
        model={1: "segmented", 2: "dual"},
        telescopes=[1, 2],
    )
    compare(command, models, tmp_path, source, telescopes=[2])


def test_wavelength_sampling_and_large_chunks(command, models, tmp_path):
    # Wavelength 0: drawn from the run's 1/lambda^2 law. Bunches of 3 photons give more output
    # rows than one chunk holds (2^16).
    photons = beam(60000, (0.0, 0.0), radius_m=6.0)
    photons["w"] = np.full(60000, 3.0, np.float32)
    photons["wavelength"] = np.zeros(60000, np.float32)
    single = dict(RUN, wavelength_min_nm=400.0, wavelength_max_nm=400.0)
    for run in (RUN, single):
        source = write_stream(tmp_path / "in.pst", {1: photons}, run=run)
        compare(command, models, tmp_path, source, threads=3)


@pytest.mark.parametrize(
    ("status", "pointing", "frame"),
    [
        (
            {"factors": "telescope_transmission", "telescope_transmission.state": "applied"},
            (0.0, 0.0),
            1,
        ),
        (None, None, 1),
        (None, (0.0, 0.0), 3),
        (
            {"factors": "mirror_reflectivity", "mirror_reflectivity.state": "envelope"},
            (0.0, 0.0),
            1,
        ),
    ],
)
def test_same_errors(command, models, tmp_path, status, pointing, frame):
    source = write_stream(
        tmp_path / "in.pst",
        {1: beam(20, (0.0, 0.0))},
        frame=frame,
        run_efficiency_statuses=[status] if status else (),
    )
    error = compare(command, models, tmp_path, source, pointing=pointing)
    assert isinstance(error, str)


def test_model_hash_is_checked(command, models, tmp_path):
    text = models["segmented"][1].read_text().replace('"focal_length_m"', '"focal_length_m" ', 1)
    changed = text.replace("0.9", "0.8", 1)
    path = tmp_path / "changed.json"
    path.write_text(changed)
    result = run_command(command, [str(path)], tmp_path / "none.pst", tmp_path / "out.pst")
    assert result.returncode == 1
    assert "does not match its sha256" in result.stderr


def test_data_after_stream_end_is_rejected(command, models, tmp_path):
    source = write_stream(tmp_path / "in.pst", {1: beam(20, (0.0, 0.0))})
    with open(source, "ab") as stream:
        stream.write(b"\x00trailing")
    with pytest.raises(ModelError, match="data after stream_end"):
        trace_stream(
            models["segmented"][0], source, tmp_path / "python.pst", seed=SEED, pointing=(0.0, 0.0)
        )
    result = run_command(command, [str(models["segmented"][1])], source, tmp_path / "native.pst")
    assert result.returncode == 1
    assert "data after stream_end" in result.stderr


def test_command_takes_a_model_per_telescope(command, models, tmp_path):
    source = write_stream(tmp_path / "in.pst", {1: beam(1400, (0.0, 0.0), radius_m=6.0)})
    output = tmp_path / "array.pst"
    argv = command_argv(
        command,
        [f"1={models['segmented'][1]}", f"2={models['dual'][1]}"],
        source,
        output,
        telescopes=None,
    )
    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "1400 photons in" in result.stderr


@pytest.mark.parametrize(
    ("models_given", "message"),
    [
        (["{path}", "2={path}"], "give one optical model for every telescope"),
        (["1={path}", "1={path}"], "more than one optical model"),
        (["0={path}"], "outside 1..65535"),
        (["LST={path}"], "--model takes PATH or TELESCOPE_ID=PATH"),
    ],
)
def test_command_rejects_inconsistent_models(command, models, tmp_path, models_given, message):
    path = models["segmented"][1]
    source = write_stream(tmp_path / "bad.in", {1: beam(10, (0.0, 0.0))})
    argv = [command, "ray-tracing", "trace-stream", "--seed", "3", "-i", str(source)]
    argv += ["-o", str(tmp_path / "bad.out")]
    for item in models_given:
        argv += ["--model", item.format(path=path)]
    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert message in result.stderr


def test_command_takes_the_telescope_transmission_option(command, models, tmp_path, obdeect_stream):
    source = write_stream(tmp_path / "opt.in", {1: beam(50, (0.0, 0.0))})
    for part in ("full", "angular"):
        output = tmp_path / f"opt.{part}.pst"
        assert (
            run_command(
                command, [str(models["segmented"][1])], source, output, transmission=part
            ).returncode
            == 0
        )
    bad = run_command(
        command, [str(models["segmented"][1])], source, tmp_path / "bad.pst", transmission="partial"
    )
    assert bad.returncode == 2
    assert "full or angular" in bad.stderr


def test_command_needs_a_pointing_for_ground_photons(command, models, tmp_path):
    source = write_stream(tmp_path / "nop.in", {1: beam(10, (0.0, 0.0))})
    result = run_command(
        command, [str(models["segmented"][1])], source, tmp_path / "nop.out", pointing=None
    )
    assert result.returncode == 1
    assert "--pointing-zenith-deg" in result.stderr


def test_command_rejects_half_a_pointing(command, models, tmp_path):
    source = write_stream(tmp_path / "half.in", {1: beam(10, (0.0, 0.0))})
    argv = [command, "ray-tracing", "trace-stream", "--model", str(models["segmented"][1])]
    argv += ["--seed", "3", "--pointing-zenith-deg", "5", "-i", str(source)]
    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert result.returncode == 1
    assert "pointing" in result.stderr
