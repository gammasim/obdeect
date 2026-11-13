"""ps2grisu: Python API and command line."""

import shutil
import subprocess
import sys
from pathlib import Path

import obdeect.stream as obdeect_stream
import numpy as np
import pytest


def _executable(name):
    local = Path(sys.executable).parent / name
    return str(local) if local.exists() else shutil.which(name)


def _write_stream(path, weight=1.0, atmosphere="applied", wavelength=400.0):
    efficiency_status = {
        "factors": "atmosphere",
        "atmosphere.state": atmosphere,
        "atmosphere.source": "pytest",
    }
    columns = {
        "x": np.array([1.5], np.float32),
        "y": np.array([0.5], np.float32),
        "cx": np.array([0.0], np.float32),
        "cy": np.array([0.5], np.float32),
        "t": np.array([12.5], np.float32),
        "w": np.array([weight], np.float32),
        "zem": np.array([20000.0], np.float32),
        "wavelength": np.array([wavelength], np.float32),
    }
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write(
            "run",
            {
                "run_number": 5,
                "observation_level_m": 2200.0,
                "telescope_x_m": np.array([10.0]),
                "telescope_y_m": np.array([20.0]),
            },
        )
        writer.write("efficiency_status", efficiency_status)
        writer.write(
            "event",
            {
                "event_number": 1,
                "primary_id": 14,
                "energy_gev": 1000.0,
                "zenith_rad": np.pi / 2,
                "azimuth_rad": 0.0,
                "first_interaction_m": 25000.0,
                "core_x_m": np.array([100.0]),
                "core_y_m": np.array([50.0]),
            },
        )
        writer.write("telescope", {"telescope_id": 1, "frame": "ground", "array_index": 0})
        writer.write_photons(1, "ground", 0, columns)
        writer.write("telescope_end", {"telescope_id": 1})
        writer.write("event_end", {})
        writer.write("run_end", {})


def _body(text):
    return text.split("* DATAF  <-- end of header flag\n", 1)[1]


def test_stream_to_grisu_text(tmp_path):
    stream = tmp_path / "run.ps"
    _write_stream(stream)
    config = tmp_path / "TLLOC.dat"
    config.write_text("* a comment\n* TLLOC 4 -20.1 -10.0 0 0 0 0\n")
    out = tmp_path / "run.grisu"
    obdeect_stream.stream_to_grisu(stream, out, config, qe=0.5, event_details=True)
    text = out.read_text()
    assert text.startswith("* HEADF  <-- Start of header flag\n")
    assert _body(text) == (
        "R 0.5\nH 2200\n"
        "S 1 -50 -100 0 -1 25000 -1 -1 -1\n"
        "C 25000 -1 1\n"
        "P -20.5 -11.5 -0.5 +0 +20000 +12.5 400 3 4\n"
    )


def test_unmatched_telescope_warns(tmp_path):
    stream = tmp_path / "run.ps"
    _write_stream(stream)
    config = tmp_path / "TLLOC.dat"
    config.write_text("* TLLOC 4 100 100 0 0 0 0\n")
    with pytest.warns(UserWarning, match="skipped 1 telescope"):
        obdeect_stream.stream_to_grisu(stream, tmp_path / "run.grisu", config)
    assert "P " not in (tmp_path / "run.grisu").read_text()


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"weight": 2.0}, "above 1"),
        ({"weight": 0.5}, "below 1"),
        ({"atmosphere": "pending"}, "atmosphere"),
        ({"wavelength": 0.0}, "wavelength"),
    ],
)
def test_errors(tmp_path, options, message):
    stream = tmp_path / "run.ps"
    _write_stream(stream, **options)
    with pytest.raises(RuntimeError, match=message):
        obdeect_stream.stream_to_grisu(stream, tmp_path / "run.grisu")


def test_expand_bunches(tmp_path):
    stream = tmp_path / "run.ps"
    _write_stream(stream, weight=3.0)
    obdeect_stream.stream_to_grisu(stream, tmp_path / "run.grisu", expand_bunches=True)
    assert (tmp_path / "run.grisu").read_text().count("\nP ") == 3


@pytest.mark.skipif(_executable("ps2grisu") is None, reason="ps2grisu is not installed")
def test_cli_matches_python(tmp_path):
    stream = tmp_path / "run.ps"
    _write_stream(stream)
    obdeect_stream.stream_to_grisu(stream, tmp_path / "py.grisu")
    result = subprocess.run(
        [_executable("ps2grisu"), str(stream)], capture_output=True, check=True, text=True
    )
    # Headers differ only in the producer line.
    assert _body(result.stdout) == _body((tmp_path / "py.grisu").read_text())
    failure = subprocess.run(
        [_executable("ps2grisu"), "--telescope-config", str(tmp_path / "missing"), str(stream)],
        capture_output=True,
        check=False,
        text=True,
    )
    assert failure.returncode == 1
    assert "cannot open telescope config" in failure.stderr
