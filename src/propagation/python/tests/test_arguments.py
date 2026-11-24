"""Argument checks of the propagation Python API: invalid input is rejected, never truncated,
cast or ignored, and a failed run leaves no output file."""

import numpy as np
import obdeect.propagation as obdeect_propagation
import pytest
from obdeect.propagation._arrays import broadcast_call, float32_columns, float_pair, integer
from test_propagate import SITE_M, TABLE_PATH, ground_stream


@pytest.fixture
def source(tmp_path):
    path = tmp_path / "in.pst"
    path.write_bytes(ground_stream(100, seed=1))
    return path


def run(source, output, **changes):
    options = {
        "transmission": TABLE_PATH,
        "site_altitude_m": SITE_M,
        "pathlen_m": 0.0,
        "seed": 11,
    }
    options.update(changes)
    return obdeect_propagation.propagate(source, output, **options)


def test_valid_run_writes_output_without_partial_file(source, tmp_path):
    output = tmp_path / "out.pst"
    run(source, output)
    assert output.exists()
    assert not (tmp_path / "out.pst.partial").exists()


@pytest.mark.parametrize(
    ("changes", "error", "message"),
    [
        ({"seed": 1.7}, TypeError, "seed must be an integer"),
        ({"seed": "11"}, TypeError, "seed must be an integer"),
        ({"threads": 2.0}, TypeError, "threads must be an integer"),
        ({"scattering_points": 4}, ValueError, "scattering_points needs scattering"),
        ({"earth_radius_m": 6.371e6}, ValueError, "earth_radius_m needs curvature"),
        ({"wavelength_range": (300.0,)}, ValueError, "exactly two elements"),
        ({"wavelength_range": (300.0, 400.0, 500.0)}, ValueError, "exactly two elements"),
        ({"envelope": "ab"}, ValueError, "envelope must be a mapping or a sequence"),
        ({"envelope": ["qe.dat"]}, ValueError, r"\(name, file\) pair"),
        ({"pathlen_m": None}, ValueError, "pathlen_m"),
    ],
)
def test_invalid_arguments_are_rejected(source, tmp_path, changes, error, message):
    output = tmp_path / "out.pst"
    with pytest.raises(error, match=message):
        run(source, output, **changes)
    assert not output.exists()


def test_failed_run_leaves_no_output(source, tmp_path):
    output = tmp_path / "out.pst"
    # The table's observation level does not match: the core rejects the run while it streams.
    with pytest.raises((ValueError, RuntimeError)):
        run(source, output, site_altitude_m=SITE_M + 1000.0)
    assert not output.exists()
    assert not (tmp_path / "out.pst.partial").exists()


def test_integer_rejects_floats_and_keeps_integers():
    assert integer(np.int64(3), "n") == 3
    with pytest.raises(TypeError, match="n must be an integer"):
        integer(3.0, "n")


def test_float_pair_needs_two_numbers():
    assert float_pair([1, 2], "r") == (1.0, 2.0)
    with pytest.raises(TypeError, match="pair of numbers"):
        float_pair("12", "r")


def test_float32_columns_reject_other_dtypes():
    columns = float32_columns({"x": np.zeros(3, dtype=np.float32)})
    assert columns["x"].dtype == np.float32
    with pytest.raises(TypeError, match="'x' must be float32"):
        float32_columns({"x": np.zeros(3)})


def test_broadcast_call_returns_array_or_tuple():
    def one(a, b):
        return a + b

    def two(a, b):
        return a + b, a - b

    assert broadcast_call(one, [1.0, 2.0], 1.0).shape == (2,)
    first, second = broadcast_call(two, [[1.0], [2.0]], [1.0, 3.0])
    assert first.shape == second.shape == (2, 2)
