"""Upward photons (cz > 0) in the atmosphere stage: loss categories, determinism, round trip."""

import tempfile
from pathlib import Path

import numpy as np
import obdeect.propagation as obdeect_propagation
import obdeect.stream as obdeect_stream
import pytest

DATA = Path(__file__).resolve().parents[2] / "cpp" / "tests" / "data"
TABLE_PATH = DATA / "atm_trans_2000_1_10_0_0_2000.dat"
PROFILE_PATH = DATA / "atmprof_ecmwf_north_winter_fixed.dat"
SITE_M = 2000.0
PLANE_Z_M = 50.0


class StreamWriter:
    """Writes through obdeect.stream.Writer into a temporary file; getvalue() returns its bytes."""

    def __init__(self):
        self._directory = tempfile.TemporaryDirectory()
        self._path = Path(self._directory.name) / "stream.pst"
        self._writer = obdeect_stream.Writer(str(self._path), "test")

    def metadata(self, kind, values=None):
        self._writer.write(kind, values or {})

    def photons(self, telescope_id, frame, columns, first_index=0):
        arrays = {name: np.asarray(values, dtype=np.float32) for name, values in columns.items()}
        self._writer.write_photons(telescope_id, frame, first_index, arrays)

    def getvalue(self):
        self._writer.close()
        try:
            return self._path.read_bytes()
        finally:
            self._directory.cleanup()


def _run(refraction_applied=0, curvature_applied=None):
    run = {
        "run_number": 1,
        "observation_level_m": SITE_M,
        "telescope_x_m": [0.0],
        "telescope_y_m": [0.0],
        "telescope_z_m": [PLANE_Z_M],
        "refraction_applied": refraction_applied,
        "wavelength_min_nm": 250.0,
        "wavelength_max_nm": 700.0,
    }
    if curvature_applied is not None:
        run["curvature_applied"] = curvature_applied
    return run


def _groups(n=10):
    """Rows of each loss category and of the two reaching categories, n photons per group."""
    c, s = 0.3, 0.2
    cz = np.sqrt(1.0 - c * c - s * s)
    # name, z [m], cx, cy, cz
    spec = {
        "upward_reach": (20.0, c, s, cz),
        "upward_above": (80.0, c, s, cz),  # above the plane: moving away
        "upward_below_ground": (-10.0, c, s, cz),
        "downward_reach": (500.0, c, s, -cz),
        "downward_above_ground_below_plane": (20.0, c, s, -cz),  # moving away
        "downward_below_ground": (-10.0, c, s, -cz),
        "horizontal": (30.0, 1.0, 0.0, 0.0),
    }
    rows = {}
    for name, (z, cx, cy, cz_value) in spec.items():
        rows[name] = {
            "x": np.full(n, 3.0),
            "y": np.full(n, -2.0),
            "z": np.full(n, z),
            "cx": np.full(n, cx),
            "cy": np.full(n, cy),
            "cz": np.full(n, cz_value),
            "t": np.zeros(n),
            "w": np.ones(n),
        }
    return rows


def _write_mixed(path, chunk_rows=None, refraction_applied=0, curvature_applied=None):
    writer = StreamWriter()
    writer.metadata("run", _run(refraction_applied, curvature_applied))
    writer.metadata("event", {"event_number": 3})
    writer.metadata("telescope", {"telescope_id": 1, "frame": 0})
    groups = _groups()
    columns = {
        name: np.concatenate([groups[group][name] for group in groups])
        for name in ("x", "y", "z", "cx", "cy", "cz", "t", "w")
    }
    rows = len(columns["x"])
    step = chunk_rows or rows
    for first in range(0, rows, step):
        end = min(rows, first + step)
        piece = {name: values[first:end] for name, values in columns.items()}
        writer.photons(1, 0, piece, first_index=first)
    writer.metadata("telescope_end", {"telescope_id": 1, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    path.write_bytes(writer.getvalue())
    return path


def test_loss_categories_are_counted_by_reason(tmp_path):
    source = _write_mixed(tmp_path / "in.pst")
    summary = obdeect_propagation.propagate(
        source, tmp_path / "out.pst", extinction=False, site_altitude_m=SITE_M, seed=1
    )
    n = 10
    assert summary["photons_in"] == 7 * n
    assert summary["photons_out"] == 2 * n
    by_reason = summary["dropped_by_reason"]
    assert by_reason == {
        "moving_away": 2 * n,
        "horizontal": n,
        "below_ground": 2 * n,
        "turns_over": 0,
        "not_bendable": 0,
        "misses_level": 0,
    }
    assert summary["dropped"] == sum(by_reason.values()) == 5 * n


def test_chunk_size_and_threads_do_not_change_the_output(tmp_path):
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 9}
    whole = _write_mixed(tmp_path / "whole.pst")
    small = _write_mixed(tmp_path / "small.pst", chunk_rows=7)
    out_whole = tmp_path / "out_whole.pst"
    out_small = tmp_path / "out_small.pst"
    out_threads = tmp_path / "out_threads.pst"
    summary_whole = obdeect_propagation.propagate(whole, out_whole, **options)
    summary_small = obdeect_propagation.propagate(small, out_small, **options)
    obdeect_propagation.propagate(whole, out_threads, threads=4, **options)
    # The records follow the input chunks, so the bytes of chunk sizes differ; the photons do not.
    assert summary_small == summary_whole
    records_whole = list(obdeect_stream.Reader(str(out_whole)))
    records_small = list(obdeect_stream.Reader(str(out_small)))
    for name in ("x", "y", "cx", "cy", "cz", "t", "zem", "w"):
        assert np.array_equal(_column(records_small, name), _column(records_whole, name))
    assert out_threads.read_bytes() == out_whole.read_bytes()


def _column(records, name):
    return np.concatenate([
        np.asarray(r["columns"][name], dtype=float) for r in records if r["type"] == "photons"
    ])


def test_cz_column_is_kept_with_its_sign_in_frame_g(tmp_path):
    source = _write_mixed(tmp_path / "in.pst")
    out = tmp_path / "out.pst"
    obdeect_propagation.propagate(source, out, extinction=False, site_altitude_m=SITE_M, seed=1)
    records = list(obdeect_stream.Reader(str(out)))
    cz = _column(records, "cz")
    # Only the two reaching groups are written: 10 upward (cz > 0) and 10 downward (cz < 0).
    assert cz.size == 20
    assert np.sum(cz > 0.0) == 10
    assert np.sum(cz < 0.0) == 10
    assert np.allclose(cx_cy_norm(records), 1.0, atol=1e-5)


def cx_cy_norm(records):
    cx = _column(records, "cx")
    cy = _column(records, "cy")
    cz = _column(records, "cz")
    return cx * cx + cy * cy + cz * cz


def test_upward_round_trip_through_to_sky(tmp_path):
    """Photons already refracted (profile, refraction_applied = 1): to_sky inverts the stage."""
    rng = np.random.default_rng(5)
    n = 200
    cx = rng.uniform(-0.3, 0.3, n)
    cy = rng.uniform(-0.3, 0.3, n)
    cz = np.sqrt(1.0 - cx**2 - cy**2)
    z = rng.uniform(0.0, PLANE_Z_M, n)  # below the plane: upward photons reach it
    writer = StreamWriter()
    writer.metadata("run", _run(refraction_applied=1))
    writer.metadata("event", {"event_number": 4})
    writer.metadata("telescope", {"telescope_id": 1, "frame": 0})
    columns = {
        "x": rng.uniform(-2000.0, 2000.0, n),
        "y": rng.uniform(-2000.0, 2000.0, n),
        "z": z,
        "cx": cx,
        "cy": cy,
        "cz": cz,
        "t": rng.uniform(0.0, 10.0, n),
        "w": np.ones(n),
    }
    writer.photons(1, 0, columns)
    writer.metadata("telescope_end", {"telescope_id": 1, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    source = tmp_path / "in.pst"
    source.write_bytes(writer.getvalue())

    ground = tmp_path / "ground.pst"
    back = tmp_path / "back.pst"
    summary = obdeect_propagation.propagate(
        source, ground, extinction=False, site_altitude_m=SITE_M, seed=1, profile=PROFILE_PATH
    )
    assert summary["photons_out"] == n
    obdeect_propagation.to_sky(ground, back, profile=PROFILE_PATH)
    x_back = _column(list(obdeect_stream.Reader(str(back))), "x")
    assert np.allclose(x_back, columns["x"], atol=1e-2)


def test_curved_mixed_directions_run_and_count(tmp_path):
    source = _write_mixed(tmp_path / "in.pst", curvature_applied=0)
    summary = obdeect_propagation.propagate(
        source,
        tmp_path / "out.pst",
        extinction=False,
        site_altitude_m=SITE_M,
        seed=1,
        curvature="spherical",
        profile=PROFILE_PATH,
    )
    assert summary["photons_in"] == 70
    assert summary["dropped"] == sum(summary["dropped_by_reason"].values())


@pytest.mark.parametrize("threads", [1, 3])
def test_threads_keep_the_mixed_output(tmp_path, threads):
    source = _write_mixed(tmp_path / "in.pst")
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 2}
    reference = tmp_path / "reference.pst"
    obdeect_propagation.propagate(source, reference, **options)
    out = tmp_path / f"out_{threads}.pst"
    obdeect_propagation.propagate(source, out, threads=threads, **options)
    assert out.read_bytes() == reference.read_bytes()
