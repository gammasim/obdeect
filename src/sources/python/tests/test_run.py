"""Runs: written photon-stream files (read back with obdeect.stream), the run description
stored in them, seeds, telescopes, and check() against the written photons."""

import json

import numpy as np
import pytest
from obdeect.sources import run
from obdeect.sources.source import photon_count

stream = pytest.importorskip("obdeect.stream")

EMISSIONS = [
    {
        "kind": "point",
        "photons": 2000.0,
        "position_m": [0.0, 0.0, 10.0],
        "axis": [0.0, 0.0, -1.0],
        "angular": {"kind": "isotropic"},
        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
        "target": {"kind": "auto"},
    },
    {
        "kind": "plane_wave",
        "photons": 300.0,
        "bunch_size": 3.0,
        "towards_source": [0.0, 0.0, 1.0],
        "disc_radius_m": 2.0,
        "disc_z_m": 0.0,
        "start_z_m": 40.0,
        "spectrum": {"kind": "uniform", "low_nm": 300.0, "high_nm": 600.0},
        "pulse": {"kind": "tophat", "width_ns": 4.0},
    },
]

RUN = {
    "run": {
        "number": 2,
        "seed": 7,
        "observation_level_m": 2150.0,
        "frame": "telescope",
        "metadata": {"note": "test"},
    },
    "telescopes": [
        {"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 12.5},
        {"id": 5, "position_m": [30.0, 0.0, 0.0], "radius_m": 6.0},
    ],
    "events": [
        {"repeat": 2, "metadata": {"label": "a"}, "emissions": [EMISSIONS[0]]},
        {"repeat": 1, "emissions": [EMISSIONS[1]]},
    ],
}


def records(path):
    return list(stream.Reader(str(path)))


def run_record(path):
    """The run record (the reader yields a stream header record first)."""
    return next(r for r in records(path) if r["type"] == "run")


def test_write_produces_a_readable_stream_with_the_run_description(tmp_path):
    path = tmp_path / "run.pst"
    summary = run.write(RUN, path)
    recs = records(path)
    types = [r["type"] for r in recs]
    assert "run" in types and types[-1] == "run_end"
    assert types.count("event") == 3
    assert types.count("telescope") == 6
    run_meta = run_record(path)["metadata"]
    description = json.loads(run_meta["run_description"])
    assert description["run"]["seed"] == 7
    assert description["events"][0]["repeat"] == 2
    assert run_meta["run_number"] == 2
    assert run_meta["note"] == "test"
    assert run_meta["n_telescopes"] == 2
    assert summary["events"] == 3


def test_run_description_is_the_normalized_run_with_sorted_keys(tmp_path):
    path = tmp_path / "run.pst"
    run.write(RUN, path)
    text = run_record(path)["metadata"]["run_description"]
    assert text == json.dumps(json.loads(text), sort_keys=True)
    assert json.loads(text)["telescopes"][0]["radius_m"] == 12.5


def test_seed_override_is_stored_and_changes_the_photons(tmp_path):
    first = tmp_path / "a.pst"
    second = tmp_path / "b.pst"
    run.write(RUN, first)
    run.write(RUN, second, seed=8)
    assert json.loads(run_record(second)["metadata"]["run_description"])["run"]["seed"] == 8
    assert not np.array_equal(_positions(first), _positions(second))


def test_invalid_seed_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="seed"):
        run.write(RUN, tmp_path / "x.pst", seed=-1)


def test_telescopes_see_the_same_photons(tmp_path):
    path = tmp_path / "run.pst"
    run.write(RUN, path)
    blocks = [r for r in records(path) if r["type"] == "photons"]
    per_telescope = {}
    for block in blocks:
        per_telescope.setdefault(block["telescope_id"], []).append(block["columns"]["x"])
    assert set(per_telescope) == {1, 5}
    assert np.array_equal(np.concatenate(per_telescope[1]), np.concatenate(per_telescope[5]))


def test_written_weights_and_rows_match_check(tmp_path):
    path = tmp_path / "run.pst"
    summary = run.write(RUN, path)
    numbers = run.check(RUN)
    assert summary["expected"] == pytest.approx(numbers["expected_weight"], rel=1e-12)
    # Point (no target inside the fiducial sphere) and plane wave keep every candidate: the
    # rows are the photon counts of each event (event numbers 1..3), for each of 2 telescopes.
    seeds = [entry["seed"] for entry in numbers["emissions"]]
    # The auto target of the point resolves to none (the point is inside the fiducial sphere).
    resolved = dict(EMISSIONS[0], target={"kind": "none"})
    counts = [
        photon_count(resolved, seed=seeds[0], event=1),
        photon_count(resolved, seed=seeds[0], event=2),
        photon_count(EMISSIONS[1], seed=seeds[1], event=3),
    ]
    assert summary["rows"] == 2 * sum(counts)
    assert summary["weight"] == pytest.approx(numbers["expected_weight"], rel=2e-2)


def test_check_counts_equal_expected_counts_and_photon_counts():
    numbers = run.check(RUN)
    assert numbers["events"] == 3
    assert numbers["telescopes"] == 2
    entries = numbers["emissions"]
    assert [e["repeat"] for e in entries] == [2, 1]
    assert entries[0]["candidates"] == pytest.approx(2000.0)
    assert entries[1]["candidates"] == pytest.approx(100.0)  # 300 photons / bunch 3
    assert entries[1]["expected_weight"] == pytest.approx(300.0, rel=1e-12)
    # Emission seeds: emission k of the run (event order) uses source_seed(run.seed, k).
    assert entries[0]["seed"] != entries[1]["seed"]
    total = sum(e["expected_weight"] * e["repeat"] for e in entries) * 2
    assert numbers["expected_weight"] == pytest.approx(total, rel=1e-12)


def test_invalid_run_is_rejected_with_the_key_path():
    broken = json.loads(json.dumps(RUN))
    broken["events"][0]["emissions"][0]["kind"] = "cones"
    with pytest.raises(ValueError, match="events\\[0\\]\\.emissions\\[0\\]"):
        run.check(broken)


def _positions(path):
    return np.concatenate([r["columns"]["x"] for r in records(path) if r["type"] == "photons"])


def test_write_rejects_runs_over_max_rows(tmp_path):
    out = tmp_path / "big.pst"
    with pytest.raises(ValueError, match="photon rows"):
        run.write(RUN, out, max_rows=1)
    assert not out.exists()
