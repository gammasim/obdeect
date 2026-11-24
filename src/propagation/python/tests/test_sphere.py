"""Sphere representation of frame G through the Python API (docs/propagation/propagation.md)."""

import tempfile
from pathlib import Path

import numpy as np
import obdeect.propagation as obdeect_propagation
import obdeect.stream as obdeect_stream
import pytest

SITE_M = 2000.0
RADIUS_M = 12.0
TELESCOPES = ((100.0, 0.0, 0.0), (-50.0, 300.0, 0.0))  # x, y, z above the observation level


class StreamWriter:
    """Writes a stream into a temporary file; getvalue() returns its bytes."""

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


def _sky_stream(rows_per_telescope=400, seed=7, radius=RADIUS_M):
    """Sky-frame photons aimed at each telescope: the line passes it at an offset of up to 1.5 R
    (about two thirds enter the sphere), emission heights above and below the telescope."""
    rng = np.random.default_rng(seed)
    writer = StreamWriter()
    writer.metadata(
        "run",
        {
            "run_number": 1,
            "observation_level_m": SITE_M,
            "telescope_x_m": [t[0] for t in TELESCOPES],
            "telescope_y_m": [t[1] for t in TELESCOPES],
            "telescope_z_m": [t[2] for t in TELESCOPES],
            "telescope_r_m": [radius, radius],
            "refraction_applied": 0,
            "wavelength_min_nm": 250.0,
            "wavelength_max_nm": 700.0,
        },
    )
    writer.metadata("event", {"event_number": 5})
    for telescope_id, (tx, ty, tz) in enumerate(TELESCOPES, start=1):
        writer.metadata("telescope", {"telescope_id": telescope_id, "frame": 0})
        n = rows_per_telescope
        cz_mag = rng.uniform(0.2, 0.999, n)
        azimuth = rng.uniform(0.0, 2.0 * np.pi, n)
        cx = np.sqrt(1.0 - cz_mag**2) * np.cos(azimuth)
        cy = np.sqrt(1.0 - cz_mag**2) * np.sin(azimuth)
        cz = rng.choice([-1.0, 1.0], n) * cz_mag
        offset = rng.uniform(-1.5 * radius, 1.5 * radius, (n, 2))
        z_emit = rng.uniform(5.0, 3000.0, n)
        # The line through Q = C + offset (at the telescope altitude) reaches z = z_emit at t.
        t_line = (z_emit - tz) / cz
        x = tx + offset[:, 0] + t_line * cx
        y = ty + offset[:, 1] + t_line * cy
        columns = {
            "x": x,
            "y": y,
            "z": z_emit,
            "cx": cx,
            "cy": cy,
            "cz": cz,
            "t": rng.uniform(0.0, 10.0, n),
            "w": np.ones(n),
        }
        writer.photons(telescope_id, 0, columns, first_index=0)
        writer.metadata("telescope_end", {"telescope_id": telescope_id, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    return writer.getvalue()


def _run(tmp_path, data, **options):
    source = tmp_path / "in.pst"
    source.write_bytes(data)
    output = tmp_path / "out.pst"
    summary = obdeect_propagation.propagate(
        source,
        output,
        extinction=False,
        site_altitude_m=SITE_M,
        seed=1,
        **options,
    )
    return summary, output


def _records(path):
    return list(obdeect_stream.Reader(str(path)))


def test_sphere_rows_entry_and_counts(tmp_path):
    summary, output = _run(tmp_path, _sky_stream(), ground_representation="sphere")
    assert summary["photons_in"] == 800
    assert summary["photons_out"] + summary["dropped"] == summary["photons_in"]
    by_reason = summary["dropped_by_reason"]
    assert "misses_sphere" in by_reason
    assert sum(by_reason.values()) == summary["dropped"]
    assert by_reason["misses_sphere"] > 0
    assert summary["photons_out"] > 0


def test_run_record_and_columns(tmp_path):
    _, output = _run(tmp_path, _sky_stream(), ground_representation="sphere")
    records = _records(output)
    run = next(r for r in records if r["type"] == "run")
    assert run["metadata"]["ground_representation"] == "sphere"
    stream = next(r for r in records if r["type"] == "stream")
    assert (stream["metadata"]["format_major"], stream["metadata"]["format_minor"]) == (1, 2)
    chunks = [r for r in records if r["type"] == "photons"]
    assert chunks
    for chunk in chunks:
        assert chunk["frame"] == 1
        columns = chunk["columns"]
        for name in ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "path"):
            assert name in columns
        assert np.all(columns["path"] >= 0.0)
        # Entry rows lie on the sphere (|r| = R), emission-inside rows inside it.
        radius = np.sqrt(columns["x"] ** 2 + columns["y"] ** 2 + columns["z"] ** 2)
        on_sphere = columns["path"] > 0.0
        assert np.all(np.abs(radius[on_sphere] - RADIUS_M) < 1e-3)
        assert np.all(radius[~on_sphere] <= RADIUS_M + 1e-3)


def test_plane_mode_is_the_default(tmp_path):
    summary, output = _run(tmp_path, _sky_stream())
    run = next(r for r in _records(output) if r["type"] == "run")
    assert "ground_representation" not in run["metadata"]
    assert summary["dropped_by_reason"]["misses_sphere"] == 0


def _horizontal_stream(rows_per_telescope=200, seed=11):
    """Horizontal sky-frame primaries (cz = 0) at 5 km above the observation level, starting 4 km
    from each telescope in random azimuths: they are scattered into the spheres (sphere runs)."""
    rng = np.random.default_rng(seed)
    writer = StreamWriter()
    writer.metadata(
        "run",
        {
            "run_number": 1,
            "observation_level_m": SITE_M,
            "telescope_x_m": [t[0] for t in TELESCOPES],
            "telescope_y_m": [t[1] for t in TELESCOPES],
            "telescope_z_m": [t[2] for t in TELESCOPES],
            "telescope_r_m": [RADIUS_M, RADIUS_M],
            "refraction_applied": 0,
            "wavelength_min_nm": 345.0,  # inside the scattering table (340 to 400 nm)
            "wavelength_max_nm": 395.0,
        },
    )
    writer.metadata("event", {"event_number": 5})
    for telescope_id, (tx, ty, _) in enumerate(TELESCOPES, start=1):
        writer.metadata("telescope", {"telescope_id": telescope_id, "frame": 0})
        n = rows_per_telescope
        azimuth = rng.uniform(0.0, 2.0 * np.pi, n)
        cx, cy = np.cos(azimuth), np.sin(azimuth)
        columns = {
            "x": tx - 4000.0 * cx + rng.uniform(-300.0, 300.0, n),
            "y": ty - 4000.0 * cy + rng.uniform(-300.0, 300.0, n),
            "z": np.full(n, 5000.0),
            "cx": cx,
            "cy": cy,
            "cz": np.zeros(n),
            "t": np.zeros(n),
            "w": np.ones(n),
        }
        writer.photons(telescope_id, 0, columns, first_index=0)
        writer.metadata("telescope_end", {"telescope_id": telescope_id, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    return writer.getvalue()


TRANSMISSION = (
    Path(__file__).resolve().parents[2]
    / "cpp"
    / "tests"
    / "data"
    / "atm_trans_2000_1_10_0_0_2000.dat"
)
SCATTERING = (
    Path(__file__).resolve().parents[2] / "cpp" / "tests" / "data" / "scattering_exponential.dat"
)


def _photon_arrays(path):
    """All photon columns of the output, in record order (for comparisons)."""
    arrays = {}
    for record in _records(path):
        if record["type"] != "photons":
            continue
        for name, values in record["columns"].items():
            arrays.setdefault(name, []).append(np.asarray(values))
    return {name: np.concatenate(parts) for name, parts in arrays.items() if parts}


def test_extinction_of_sphere_rows_is_applied(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(_sky_stream())
    summary = {}
    outputs = {}
    for threads in (1, 3):
        output = tmp_path / f"out{threads}.pst"
        summary[threads] = obdeect_propagation.propagate(
            source,
            output,
            pathlen_m=0.0,
            transmission=TRANSMISSION,
            site_altitude_m=SITE_M,
            seed=1,
            ground_representation="sphere",
            threads=threads,
        )
        outputs[threads] = output
    assert summary[1]["photons_out"] > 0
    assert summary[1]["photons_out"] == summary[3]["photons_out"]
    efficiency = [r for r in _records(outputs[1]) if r["type"] == "efficiency_status"]
    assert efficiency
    assert efficiency[-1]["metadata"]["factors"] == "atmosphere,wavelength"
    assert efficiency[-1]["metadata"]["atmosphere.state"] == "applied"
    # Thinned output: the photons of both thread counts are identical.
    one, three = _photon_arrays(outputs[1]), _photon_arrays(outputs[3])
    assert set(one) == set(three)
    for name in one:
        np.testing.assert_array_equal(one[name], three[name])
    # Weighted output keeps the extinction in the weights (the atmosphere factor is below 1).
    weighted = tmp_path / "weighted.pst"
    obdeect_propagation.propagate(
        source,
        weighted,
        pathlen_m=0.0,
        transmission=TRANSMISSION,
        site_altitude_m=SITE_M,
        seed=1,
        ground_representation="sphere",
        mode="weighted",
    )
    w = _photon_arrays(weighted)["w"]
    assert np.all(w >= 0.0)
    assert np.any(w < 1.0)
    assert np.all(w <= 1.0 + 1e-6)


def test_horizontal_primaries_scatter_in_sphere_runs(tmp_path):
    source = tmp_path / "horizontal.pst"
    source.write_bytes(_horizontal_stream())
    output = tmp_path / "out.pst"
    # Weighted output: the scattered weights are tiny, thinning would drop them all.
    summary = obdeect_propagation.propagate(
        source,
        output,
        pathlen_m=0.0,
        transmission=TRANSMISSION,
        scattering=SCATTERING,
        scattering_points=4,
        site_altitude_m=SITE_M,
        seed=3,
        ground_representation="sphere",
        mode="weighted",
    )
    assert summary["photons_out"] > 0
    efficiency = [r for r in _records(output) if r["type"] == "efficiency_status"]
    assert "scattering" in efficiency[-1]["metadata"]["factors"]
    arrays = _photon_arrays(output)
    # Sphere rows of the scattered photons: on or inside the sphere (path >= 0), direction unit.
    radius = np.sqrt(arrays["x"] ** 2 + arrays["y"] ** 2 + arrays["z"] ** 2)
    assert np.all(radius <= RADIUS_M + 1e-3)
    np.testing.assert_allclose(
        arrays["cx"] ** 2 + arrays["cy"] ** 2 + arrays["cz"] ** 2, 1.0, atol=1e-5
    )
    # The same run with 2 threads gives the same photons.
    output2 = tmp_path / "out2.pst"
    obdeect_propagation.propagate(
        source,
        output2,
        pathlen_m=0.0,
        transmission=TRANSMISSION,
        scattering=SCATTERING,
        scattering_points=4,
        site_altitude_m=SITE_M,
        seed=3,
        ground_representation="sphere",
        mode="weighted",
        threads=2,
    )
    other = _photon_arrays(output2)
    assert set(other) == set(arrays)
    for name in arrays:
        np.testing.assert_array_equal(arrays[name], other[name])


def test_invalid_ground_representation_and_missing_radius(tmp_path):
    with pytest.raises(ValueError, match="ground representation must be plane or sphere"):
        _run(tmp_path, _sky_stream(), ground_representation="cube")
    stream = _sky_stream()
    source = tmp_path / "in.pst"
    source.write_bytes(stream)
    # A run without telescope_r_m is rejected in sphere mode.
    writer = StreamWriter()
    writer.metadata(
        "run",
        {
            "run_number": 1,
            "observation_level_m": SITE_M,
            "telescope_x_m": [100.0],
            "telescope_y_m": [0.0],
            "telescope_z_m": [0.0],
            "refraction_applied": 0,
            "wavelength_min_nm": 250.0,
            "wavelength_max_nm": 700.0,
        },
    )
    writer.metadata("run_end")
    no_radius = tmp_path / "no_radius.pst"
    no_radius.write_bytes(writer.getvalue())
    with pytest.raises(ValueError, match="telescope_r_m"):
        obdeect_propagation.propagate(
            no_radius,
            tmp_path / "out.pst",
            extinction=False,
            site_altitude_m=SITE_M,
            seed=1,
            ground_representation="sphere",
        )
