"""CORSIKA <-> photon stream <-> sim_telarray conversions through the Python API and the CLIs."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import obdeect.stream as obdeect_stream
import numpy as np
import pytest

COLUMNS = ("x", "y", "cx", "cy", "t", "w", "zem", "wavelength")


def _executable(name):
    """Installed tool next to the interpreter (a venv's bin) or on PATH; None if absent."""
    local = Path(sys.executable).parent / name
    if local.exists():
        return str(local)
    return shutil.which(name)


def _photons(rows, telescope_offset=0.0):
    """Photon columns with exact cm values (multiples of 0.25 m) and the sim_telarray fields."""
    index = np.arange(rows, dtype=np.float32)
    return {
        "x": (np.float32(0.25) * (index % 40) - np.float32(5.0)).astype(np.float32),
        "y": (np.float32(0.25) * (index % 24) + np.float32(telescope_offset)).astype(np.float32),
        "cx": (np.float32(0.01) * (index % 7)).astype(np.float32),
        "cy": (np.float32(-0.02) * (index % 5)).astype(np.float32),
        "t": index * np.float32(0.5),
        "w": np.ones(rows, dtype=np.float32),
        "zem": np.full(rows, 20000.0, dtype=np.float32),
        "wavelength": np.full(rows, 400.0, dtype=np.float32),
    }


def _write_light_source_stream(
    path, rows=50, run_efficiency_status=None, telescope_efficiency_status=None
):
    """A stream without CORSIKA blobs: telescope 1 ground frame, telescope 2 camera frame.

    ``run_efficiency_status`` is written before the event; ``telescope_efficiency_status``
    follows telescope 2's record.
    """
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write(
            "run",
            {
                "run_number": 4,
                "observation_level_m": 2200.0,
                "telescope_x_m": np.array([0.0, 100.0]),
                "telescope_y_m": np.array([0.0, 0.0]),
                "telescope_z_m": np.array([0.0, 0.0]),
                "telescope_r_m": np.array([12.0, 12.0]),
            },
        )
        if run_efficiency_status is not None:
            writer.write("efficiency_status", run_efficiency_status)
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": "ground", "array_index": 0})
        writer.write_photons(1, "ground", 0, _photons(rows))
        writer.write("telescope_end", {"telescope_id": 1, "photons": float(rows)})
        writer.write("telescope", {"telescope_id": 2, "frame": "camera", "array_index": 0})
        if telescope_efficiency_status is not None:
            writer.write("efficiency_status", telescope_efficiency_status)
        writer.write_photons(2, "camera", 0, _photons(rows, telescope_offset=1.0))
        writer.write("telescope_end", {"telescope_id": 2, "photons": float(rows)})
        writer.write("event_end", {})
        writer.write("run_end", {})


def _photons_by_telescope(path):
    """Photon columns per telescope id from a stream file."""
    result = {}
    for record in obdeect_stream.Reader(str(path)):
        if record["type"] == "photons":
            result.setdefault(record["telescope_id"], []).append(record)
    return result


def _column(records, name):
    return np.concatenate([record["columns"][name] for record in records])


def test_stream_to_simtel_and_back(tmp_path):
    stream = tmp_path / "light.ps"
    simtel = tmp_path / "light.simtel"
    back = tmp_path / "back.ps"
    _write_light_source_stream(stream)
    obdeect_stream.stream_to_simtel(stream, simtel, {2: (1, 0.0)})
    assert simtel.read_bytes()[:4] == bytes.fromhex("378a1fd4")  # EventIO sync tag
    obdeect_stream.corsika_to_stream(simtel, back)

    by_telescope = _photons_by_telescope(back)
    assert set(by_telescope) == {1, 2}
    ground = _photons_by_telescope(stream)[1]
    np.testing.assert_array_equal(_column(by_telescope[1], "x"), _column(ground, "x"))
    np.testing.assert_array_equal(_column(by_telescope[1], "y"), _column(ground, "y"))
    np.testing.assert_array_equal(_column(by_telescope[1], "w"), _column(ground, "w"))
    assert len(_column(by_telescope[2], "x")) == 50


@pytest.mark.parametrize("camera_map", ["dict", "file"])
def test_camera_frame_conversion(tmp_path, camera_map):
    stream = tmp_path / "light.ps"
    simtel = tmp_path / "light.simtel"
    back = tmp_path / "back.ps"
    _write_light_source_stream(stream)
    original = _photons_by_telescope(stream)[2]
    if camera_map == "dict":
        obdeect_stream.stream_to_simtel(stream, simtel, {2: (1, 0.0)})
    else:
        map_file = tmp_path / "camera.txt"
        map_file.write_text("# telescope mirror_class rotation\n2 1 0.0\n")
        obdeect_stream.stream_to_simtel(stream, simtel, map_file)
    obdeect_stream.corsika_to_stream(simtel, back)
    converted = _photons_by_telescope(back)[2]
    # Single mirror: (x, y) -> (x, -y) and (cx, cy) -> (cx, -cy), in metres.
    np.testing.assert_allclose(_column(converted, "x"), _column(original, "x"), atol=1e-6)
    np.testing.assert_allclose(_column(converted, "y"), -_column(original, "y"), atol=1e-6)
    np.testing.assert_allclose(_column(converted, "cy"), -_column(original, "cy"), atol=1e-6)


def _write_flasher_stream(path):
    """One telescope-frame photon: emission point (0.1, 0.2, 0.5) m, direction (0.6, 0, -0.8)."""

    def column(value):
        return np.array([value], dtype=np.float32)

    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write(
            "run",
            {
                "run_number": 4,
                "observation_level_m": 2200.0,
                "telescope_x_m": np.array([0.0]),
                "telescope_y_m": np.array([0.0]),
                "telescope_z_m": np.array([0.0]),
                "telescope_r_m": np.array([1.2]),
            },
        )
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": "telescope", "array_index": 0})
        writer.write_photons(
            1,
            "telescope",
            0,
            {
                "x": column(0.1),
                "y": column(0.2),
                "z": column(0.5),
                "cx": column(0.6),
                "cy": column(0.0),
                "cz": column(-0.8),
                "t": column(3.0),
                "w": column(1.0),
            },
        )
        writer.write("telescope_end", {"telescope_id": 1, "photons": 1.0})
        writer.write("event_end", {})
        writer.write("run_end", {})


def test_telescope_frame_conversion(tmp_path):
    stream = tmp_path / "flasher.ps"
    simtel = tmp_path / "flasher.simtel"
    back = tmp_path / "back.ps"
    _write_flasher_stream(stream)
    obdeect_stream.stream_to_simtel(stream, simtel, {1: (2, 0.0)})
    obdeect_stream.corsika_to_stream(simtel, back)
    converted = _photons_by_telescope(back)[1]
    # Crossing of the plane z = 0 after 0.625 m: x = 0.1 + 0.625 * 0.6.
    np.testing.assert_allclose(_column(converted, "x"), [0.475], atol=1e-6)
    np.testing.assert_allclose(_column(converted, "y"), [0.2], atol=1e-6)
    np.testing.assert_allclose(_column(converted, "t"), [3.0 + 0.625 / 0.299792458], atol=1e-4)
    np.testing.assert_allclose(_column(converted, "zem"), [2200.5], atol=1e-3)


def test_telescope_frame_needs_dual_mirror_entry(tmp_path):
    stream = tmp_path / "flasher.ps"
    _write_flasher_stream(stream)
    with pytest.raises(RuntimeError, match="camera map entry"):
        obdeect_stream.stream_to_simtel(stream, tmp_path / "a.simtel")
    with pytest.raises(RuntimeError, match="dual-mirror"):
        obdeect_stream.stream_to_simtel(stream, tmp_path / "b.simtel", {1: (1, 0.0)})


def test_rejects_sky_frame(tmp_path):
    stream = tmp_path / "sky.ps"
    with obdeect_stream.Writer(str(stream), "pytest") as writer:
        writer.write(
            "run",
            {
                "run_number": 1,
                "observation_level_m": 2200.0,
                "telescope_x_m": np.array([0.0]),
                "telescope_y_m": np.array([0.0]),
                "telescope_z_m": np.array([0.0]),
                "telescope_r_m": np.array([12.0]),
            },
        )
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1})
        writer.write_photons(1, "sky", 0, _photons(3))
        writer.write("telescope_end", {"telescope_id": 1, "photons": 3.0})
        writer.write("event_end", {})
        writer.write("run_end", {})
    with pytest.raises(RuntimeError, match="sky frame"):
        obdeect_stream.stream_to_simtel(stream, tmp_path / "sky.simtel")


def test_rejects_invalid_mirror_class(tmp_path):
    stream = tmp_path / "light.ps"
    _write_light_source_stream(stream)
    with pytest.raises(ValueError):
        obdeect_stream.stream_to_simtel(stream, tmp_path / "bad.simtel", {2: (3, 0.0)})


def _efficiency_status(factor, state):
    return {"factors": factor, f"{factor}.state": state}


# Efficiency status of a camera-frame telescope (docs/usage.md): accepted and rejected factors.
ACCEPTED_EFFICIENCY_STATUSES = [
    _efficiency_status("qe", "pending"),
    _efficiency_status("atmosphere", "applied"),
    _efficiency_status("telescope_transmission_constant", "applied"),
]
REJECTED_EFFICIENCY_STATUSES = [
    _efficiency_status("telescope_transmission_constant", "pending"),
    _efficiency_status("atmosphere", "pending"),
    _efficiency_status("qe", "applied"),
    _efficiency_status("qe", "envelope"),
    _efficiency_status("made_up", "applied"),
]


@pytest.mark.parametrize("efficiency_status", ACCEPTED_EFFICIENCY_STATUSES)
def test_efficiency_status_accepted_in_camera_frame(tmp_path, efficiency_status):
    stream = tmp_path / "light.ps"
    simtel = tmp_path / "light.simtel"
    _write_light_source_stream(stream, telescope_efficiency_status=efficiency_status)
    obdeect_stream.stream_to_simtel(stream, simtel, {2: (1, 0.0)})
    assert simtel.stat().st_size > 0


@pytest.mark.parametrize("efficiency_status", REJECTED_EFFICIENCY_STATUSES)
def test_efficiency_status_rejected_in_camera_frame(tmp_path, efficiency_status):
    stream = tmp_path / "light.ps"
    _write_light_source_stream(stream, telescope_efficiency_status=efficiency_status)
    with pytest.raises(RuntimeError, match=efficiency_status["factors"]):
        obdeect_stream.stream_to_simtel(stream, tmp_path / "bad.simtel", {2: (1, 0.0)})


def test_efficiency_status_pending_accepted_with_allow_pending(tmp_path):
    stream = tmp_path / "light.ps"
    simtel = tmp_path / "light.simtel"
    constant = _efficiency_status("telescope_transmission_constant", "pending")
    _write_light_source_stream(stream, telescope_efficiency_status=constant)
    obdeect_stream.stream_to_simtel(
        stream, simtel, {2: (1, 0.0)}, allow_pending=["telescope_transmission_constant"]
    )
    assert simtel.stat().st_size > 0
    with pytest.raises(TypeError):
        obdeect_stream.stream_to_simtel(stream, simtel, {2: (1, 0.0)}, allow_pending="qe")


def test_efficiency_status_scaled_factor_accepts_envelope(tmp_path):
    # A QE table divided by the QE envelope: sim_telarray applies f/f_env.
    stream = tmp_path / "light.ps"
    simtel = tmp_path / "light.simtel"
    _write_light_source_stream(
        stream,
        run_efficiency_status=_efficiency_status("qe", "envelope"),
        telescope_efficiency_status=_efficiency_status("qe", "envelope"),
    )
    with pytest.raises(RuntimeError, match="--scaled-factor qe"):
        obdeect_stream.stream_to_simtel(stream, tmp_path / "bad.simtel", {2: (1, 0.0)})
    obdeect_stream.stream_to_simtel(stream, simtel, {2: (1, 0.0)}, scaled_factors=["qe"])
    assert simtel.stat().st_size > 0
    with pytest.raises(TypeError, match="scaled_factors"):
        obdeect_stream.stream_to_simtel(stream, simtel, {2: (1, 0.0)}, scaled_factors="qe")


@pytest.mark.parametrize(
    "efficiency_status",
    [
        _efficiency_status("qe", "pending"),
        _efficiency_status("qe", "applied"),
        _efficiency_status("atmosphere", "applied"),
    ],
)
def test_efficiency_status_scaled_factor_needs_envelope_applied_by_sim_telarray(
    tmp_path, efficiency_status
):
    stream = tmp_path / "light.ps"
    _write_light_source_stream(
        stream,
        run_efficiency_status=_efficiency_status(efficiency_status["factors"], "envelope"),
        telescope_efficiency_status=efficiency_status,
    )
    name = efficiency_status["factors"]
    with pytest.raises(RuntimeError, match="given as scaled"):
        obdeect_stream.stream_to_simtel(
            stream, tmp_path / "bad.simtel", {2: (1, 0.0)}, scaled_factors=[name]
        )


def test_run_efficiency_status_applies_to_camera_telescope(tmp_path):
    # A run-level factor applies to every telescope; mirror reflectivity is applied by
    # sim_telarray in ground frame only, so the camera telescope must reject it.
    stream = tmp_path / "light.ps"
    _write_light_source_stream(
        stream, run_efficiency_status=_efficiency_status("mirror_reflectivity", "pending")
    )
    with pytest.raises(RuntimeError, match="mirror_reflectivity"):
        obdeect_stream.stream_to_simtel(stream, tmp_path / "bad.simtel", {2: (1, 0.0)})


def test_rejects_truncated_stream(tmp_path):
    stream = tmp_path / "light.ps"
    _write_light_source_stream(stream)
    truncated = tmp_path / "truncated.ps"
    truncated.write_bytes(stream.read_bytes()[:-40])
    with pytest.raises(RuntimeError):
        obdeect_stream.stream_to_simtel(truncated, tmp_path / "truncated.simtel")


@pytest.mark.skipif(
    _executable("ps2simtel") is None or _executable("corsika2ps") is None,
    reason="corsika2ps and ps2simtel are not installed",
)
def test_cli_matches_python_api(tmp_path):
    stream = tmp_path / "light.ps"
    camera = tmp_path / "camera.txt"
    api_out = tmp_path / "api.simtel"
    cli_out = tmp_path / "cli.simtel"
    _write_light_source_stream(stream)
    camera.write_text("2 1 0.0\n")
    obdeect_stream.stream_to_simtel(stream, api_out, {2: (1, 0.0)})
    subprocess.run(
        [_executable("ps2simtel"), "--camera-map", str(camera), "-o", str(cli_out), str(stream)],
        check=True,
        capture_output=True,
    )
    assert cli_out.read_bytes() == api_out.read_bytes()

    # corsika2ps reads the file from stdin and writes the stream to stdout; pipe it back.
    with cli_out.open("rb") as source:
        piped = subprocess.run(
            [_executable("corsika2ps")], stdin=source, check=True, capture_output=True
        ).stdout
    assert piped[:4] == b"GSTR"
    back = tmp_path / "piped.ps"
    back.write_bytes(piped)
    assert set(_photons_by_telescope(back)) == {1, 2}
    assert _column(_photons_by_telescope(back)[1], "x").size == 50


@pytest.mark.skipif(
    _executable("ps2simtel") is None,
    reason="ps2simtel is not installed",
)
def test_cli_allow_pending_option(tmp_path):
    stream = tmp_path / "light.ps"
    camera = tmp_path / "camera.txt"
    camera.write_text("2 1 0.0\n")
    constant = _efficiency_status("telescope_transmission_constant", "pending")
    _write_light_source_stream(stream, telescope_efficiency_status=constant)
    base = [_executable("ps2simtel"), "--camera-map", str(camera), "-o", str(tmp_path / "x")]
    failed = subprocess.run([*base, str(stream)], capture_output=True)
    assert failed.returncode != 0
    assert b"telescope_transmission_constant" in failed.stderr
    subprocess.run(
        [
            *base[:-2],
            "-o",
            str(tmp_path / "ok.simtel"),
            "--allow-pending",
            "telescope_transmission_constant",
            str(stream),
        ],
        check=True,
        capture_output=True,
    )
    assert (tmp_path / "ok.simtel").stat().st_size > 0


@pytest.mark.parametrize("camera", [{0: (1, 0.0)}, {2: (2, float("nan"))}, {2: (2, float("inf"))}])
def test_rejects_invalid_camera_mapping(tmp_path, camera):
    stream = tmp_path / "light.ps"
    _write_light_source_stream(stream)
    with pytest.raises((ValueError, RuntimeError)):
        obdeect_stream.stream_to_simtel(stream, tmp_path / "bad.simtel", camera)


@pytest.mark.parametrize("operation", ["stream_to_simtel", "corsika_to_stream"])
def test_conversion_reports_output_failure(tmp_path, operation):
    stream = tmp_path / "light.ps"
    simtel = tmp_path / "light.simtel"
    _write_light_source_stream(stream)
    obdeect_stream.stream_to_simtel(stream, simtel, {2: (1, 0.0)})
    source = stream if operation == "stream_to_simtel" else simtel
    code = """
import resource
import signal
import sys
import obdeect.stream as obdeect_stream
obdeect_stream.__path__.append(sys.argv[4])
signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
resource.setrlimit(resource.RLIMIT_FSIZE, (0, 0))
try:
    if sys.argv[1] == 'stream_to_simtel':
        obdeect_stream.stream_to_simtel(sys.argv[2], sys.argv[3], {2: (1, 0.0)})
    else:
        obdeect_stream.corsika_to_stream(sys.argv[2], sys.argv[3])
except RuntimeError:
    sys.exit(0)
sys.exit(1)
"""
    result = subprocess.run(
        [
            sys.executable,
            "-S",
            "-c",
            code,
            operation,
            str(source),
            str(tmp_path / "full"),
            str(Path(obdeect_stream._core.__file__).parent),
        ],
        capture_output=True,
        check=False,
        env={
            **os.environ,
            "PYTHONPATH": os.pathsep.join([str(Path(obdeect_stream.__file__).parent.parent), *sys.path]),
        },
    )
    assert result.returncode == 0, result.stderr.decode()


@pytest.mark.skipif(
    _executable("ps2simtel") is None,
    reason="ps2simtel is not installed",
)
def test_cli_scaled_factor_option(tmp_path):
    stream = tmp_path / "light.ps"
    camera = tmp_path / "camera.txt"
    camera.write_text("2 1 0.0\n")
    _write_light_source_stream(
        stream,
        run_efficiency_status=_efficiency_status("qe", "envelope"),
        telescope_efficiency_status=_efficiency_status("qe", "envelope"),
    )
    out = tmp_path / "ok.simtel"
    base = [_executable("ps2simtel"), "--camera-map", str(camera), "-o", str(out)]
    failed = subprocess.run([*base, str(stream)], capture_output=True)
    assert failed.returncode != 0
    assert b"--scaled-factor qe" in failed.stderr
    subprocess.run([*base, "--scaled-factor", "qe", str(stream)], check=True, capture_output=True)
    assert out.stat().st_size > 0
