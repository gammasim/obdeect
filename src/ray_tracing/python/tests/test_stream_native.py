"""The native trace-stream program writes the same stream as the Python trace_stream.

Each case traces one input with both and compares every record (columns bit for bit, metadata
exactly; only the version in the producer and provenance strings may differ), or checks that both
reject the input with the same message.
"""

import os
import re
import subprocess

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

PROGRAM = "obdeect-ray-tracing-trace-stream"


def _native_program():
    path = os.environ.get("OBDEECT_TRACE_STREAM")
    if path:
        return path
    from obdeect.stream.cli import find_program

    return find_program("obdeect.ray_tracing", PROGRAM)


@pytest.fixture(scope="module")
def native():
    path = _native_program()
    if path is None:
        pytest.skip(f"{PROGRAM} is not installed (set OBDEECT_TRACE_STREAM to a build)")
    return path


@pytest.fixture
def obdeect_stream():
    return pytest.importorskip("obdeect.stream")


@pytest.fixture(scope="module")
def models(tmp_path_factory):
    folder = tmp_path_factory.mktemp("models")
    out = {}
    for name, build in (("segmented", segmented), ("dual", dual)):
        records, files = build()
        model = OpticalModel.from_records(records, files.__getitem__, configuration_seed=1)
        path = folder / f"{name}.json"
        model.save(path)
        out[name] = (OpticalModel.load(path), path)
    return out


def _normalised(value):
    if isinstance(value, str):
        return re.sub(r"obdeect \d+\.\d+\.\d+\S*", "obdeect VERSION", value)
    return value


def _records(path, obdeect_stream):
    out = []
    for record in obdeect_stream.Reader(str(path)):
        record = dict(record)
        if "metadata" in record:
            record["metadata"] = {k: _normalised(v) for k, v in record["metadata"].items()}
        out.append(record)
    return out


def _assert_same_stream(python_out, native_out, obdeect_stream):
    left, right = _records(python_out, obdeect_stream), _records(native_out, obdeect_stream)
    assert [r["type"] for r in left] == [r["type"] for r in right]
    for a, b in zip(left, right, strict=True):
        assert a.keys() == b.keys()
        for key in a:
            if key == "columns":
                assert a[key].keys() == b[key].keys()
                for name in a[key]:
                    assert a[key][name].tobytes() == b[key][name].tobytes(), name
            elif key == "metadata":
                assert a[key].keys() == b[key].keys()
                for name, value in a[key].items():
                    assert np.array_equal(np.asarray(value), np.asarray(b[key][name])), name
            else:
                assert a[key] == b[key], key


def _native_argv(
    native, model_args, source, output, *, telescopes, pointing, transmission, threads
):
    argv = [native, "--seed", str(SEED), "--threads", str(threads)]
    for item in model_args:
        argv += ["--model", item]
    if telescopes is not None:
        argv += ["--telescopes", ",".join(map(str, telescopes))]
    if pointing is not None:
        argv += ["--pointing-zenith-deg", repr(pointing[0])]
        argv += ["--pointing-azimuth-deg", repr(pointing[1])]
    return argv + ["--telescope-transmission", transmission, "-i", str(source), "-o", str(output)]


def compare(
    native,
    models,
    tmp_path,
    obdeect_stream,
    source,
    *,
    model="segmented",
    telescopes=None,
    pointing=(0.0, 0.0),
    transmission="full",
    threads=1,
):
    """Trace `source` with both; equal streams, or the same error from both."""
    if isinstance(model, dict):
        python_model = {tel: models[name][0] for tel, name in model.items()}
        model_args = [f"{tel}={models[name][1]}" for tel, name in model.items()]
    else:
        python_model = models[model][0]
        model_args = [str(models[model][1])]
    python_out, native_out = tmp_path / "python.pst", tmp_path / "native.pst"
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
    except (ModelError, ValueError) as error:
        python_error = str(error)
    argv = _native_argv(
        native,
        model_args,
        source,
        native_out,
        telescopes=telescopes,
        pointing=pointing,
        transmission=transmission,
        threads=threads,
    )
    result = subprocess.run(argv, capture_output=True, text=True, check=False)
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
    _assert_same_stream(python_out, native_out, obdeect_stream)
    return python_counts


def test_ground_frame_on_and_off_axis(native, models, tmp_path, obdeect_stream):
    photons = beam(4000, (0.02, -0.01), radius_m=6.0)
    photons["zem"] = np.linspace(8000.0, 12000.0, 4000, dtype=np.float32)
    source = write_stream(tmp_path / "in.pst", {1: photons}, blob=True)
    counts = compare(native, models, tmp_path, obdeect_stream, source)
    assert counts["photons_out"] > 0


@pytest.mark.parametrize("cz", [-0.000614837, 0.000614837, -0.05])
def test_near_horizontal_cz_validation(native, models, tmp_path, obdeect_stream, cz):
    photons = beam(1, (0.989117742, -0.147124648))
    photons["cz"] = np.array([cz], dtype=np.float32)
    source = write_stream(tmp_path / "in.pst", {1: photons})
    result = compare(native, models, tmp_path, obdeect_stream, source)
    if 0 < abs(cz) < 0.01:
        # Nearly horizontal photons, either way up, are traced (upwards ones are not detected).
        assert result["photons_in"] == 1
    else:
        assert "cz" in result


def test_mixed_directions_match_python(native, models, tmp_path, obdeect_stream):
    down = with_direction(beam(3000, (0.0, 0.0), radius_m=0.5), 0.0, 0.0, -1.0)
    up = with_direction(beam(1500, (0.01, 0.0), radius_m=0.5, seed=8), 0.01, 0.0, +1.0)
    photons = {name: np.concatenate([down[name], up[name]]) for name in down}
    source = write_stream(tmp_path / "in.pst", {1: photons})
    counts = compare(native, models, tmp_path, obdeect_stream, source)
    assert counts["photons_in"] == 4500
    assert counts["photons_out"] > 0


def test_horizontal_photons_error_matches_python(native, models, tmp_path, obdeect_stream):
    photons = with_direction(beam(5, (1.0, 0.0)), 1.0, 0.0, 1.0)
    photons["cz"] = np.zeros(5, np.float32)
    source = write_stream(tmp_path / "in.pst", {1: photons})
    error = compare(native, models, tmp_path, obdeect_stream, source)
    assert isinstance(error, str) and "horizontal photons" in error


def test_stdin_matches_file_input(native, models, tmp_path):
    source = write_stream(tmp_path / "in.pst", {1: beam(4000, (0.0, 0.0), radius_m=6.0)})
    args = dict(telescopes=None, pointing=(0.0, 0.0), transmission="full", threads=1)
    model_args = [str(models["segmented"][1])]
    file_result = subprocess.run(
        _native_argv(native, model_args, source, "-", **args), capture_output=True, check=False
    )
    stdin_result = subprocess.run(
        _native_argv(native, model_args, "-", "-", **args),
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
def test_tilted_pointing(native, models, tmp_path, obdeect_stream, pointing):
    for frame, photons in zip(
        (2, 1), parallel_beam_in_both_frames(models["segmented"][0], *pointing, 3000, 5.0)
    ):
        source = write_stream(tmp_path / f"in{frame}.pst", {1: photons}, frame=frame)
        compare(native, models, tmp_path, obdeect_stream, source, pointing=pointing)


def test_run_pointing_and_dual_mirror(native, models, tmp_path, obdeect_stream):
    run = dict(RUN, telescope_pointing_zenith_rad=0.3, telescope_pointing_azimuth_rad=1.1)
    photons = parallel_beam_in_both_frames(
        models["dual"][0], 0.3 * 180 / np.pi, 1.1 * 180 / np.pi, 3000, 2.0
    )[1]
    source = write_stream(tmp_path / "in.pst", {1: photons}, run=run)
    compare(native, models, tmp_path, obdeect_stream, source, model="dual", pointing=None)


def test_envelope_status_and_other_factors(native, models, tmp_path, obdeect_stream):
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
    compare(native, models, tmp_path, obdeect_stream, source)


@pytest.mark.parametrize("transmission", ["full", "angular"])
def test_transmission_parts(native, models, tmp_path, obdeect_stream, transmission):
    source = write_stream(tmp_path / "in.pst", {1: beam(3000, (0.01, 0.0), radius_m=6.0)})
    compare(native, models, tmp_path, obdeect_stream, source, transmission=transmission)


def test_per_telescope_models_and_selection(native, models, tmp_path, obdeect_stream):
    photons = {
        1: beam(2000, (0.0, 0.0), radius_m=6.0),
        2: beam(2000, (0.0, 0.0), radius_m=2.0),
        3: beam(50, (0.0, 0.0)),
    }
    source = write_stream(tmp_path / "in.pst", photons)
    compare(
        native,
        models,
        tmp_path,
        obdeect_stream,
        source,
        model={1: "segmented", 2: "dual"},
        telescopes=[1, 2],
    )
    compare(native, models, tmp_path, obdeect_stream, source, telescopes=[2])


def test_wavelength_sampling_and_large_chunks(native, models, tmp_path, obdeect_stream):
    # Wavelength 0: drawn from the run's 1/lambda^2 law. Bunches of 3 photons give more output
    # rows than one chunk holds (2^16).
    photons = beam(60000, (0.0, 0.0), radius_m=6.0)
    photons["w"] = np.full(60000, 3.0, np.float32)
    photons["wavelength"] = np.zeros(60000, np.float32)
    single = dict(RUN, wavelength_min_nm=400.0, wavelength_max_nm=400.0)
    for run in (RUN, single):
        source = write_stream(tmp_path / "in.pst", {1: photons}, run=run)
        compare(native, models, tmp_path, obdeect_stream, source, threads=3)


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
def test_same_errors(native, models, tmp_path, obdeect_stream, status, pointing, frame):
    source = write_stream(
        tmp_path / "in.pst",
        {1: beam(20, (0.0, 0.0))},
        frame=frame,
        run_efficiency_statuses=[status] if status else (),
    )
    error = compare(native, models, tmp_path, obdeect_stream, source, pointing=pointing)
    assert isinstance(error, str)


def test_model_hash_is_checked(native, models, tmp_path):
    text = models["segmented"][1].read_text().replace('"focal_length_m"', '"focal_length_m" ', 1)
    changed = text.replace("0.9", "0.8", 1)
    path = tmp_path / "changed.json"
    path.write_text(changed)
    result = subprocess.run(
        [native, "--model", str(path), "--seed", "1", "-i", "-", "-o", "-"],
        capture_output=True,
        text=True,
        input="",
        check=False,
    )
    assert result.returncode == 1
    assert "does not match its sha256" in result.stderr


def test_data_after_stream_end_is_rejected(native, models, tmp_path):
    source = write_stream(tmp_path / "in.pst", {1: beam(20, (0.0, 0.0))})
    with open(source, "ab") as stream:
        stream.write(b"\x00trailing")
    with pytest.raises(Exception):
        trace_stream(models["segmented"][0], source, tmp_path / "python.pst", seed=SEED)
    out = tmp_path / "native.pst"
    argv = _native_argv(
        native,
        [str(models["segmented"][1])],
        source,
        out,
        telescopes=None,
        pointing=(0.0, 0.0),
        transmission="full",
        threads=1,
    )
    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert result.returncode == 1
    assert "data after stream_end" in result.stderr
