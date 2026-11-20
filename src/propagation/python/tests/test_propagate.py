"""Atmosphere stage through the Python API and the obdeect propagation run executable."""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import obdeect.propagation as obdeect_propagation
import obdeect.stream as obdeect_stream
import pytest
from obdeect.stream.cli import find_program

DATA = Path(__file__).resolve().parents[2] / "cpp" / "tests" / "data"
TABLE_PATH = DATA / "atm_trans_2000_1_10_0_0_2000.dat"
SITE_M = 2000.0


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


def read_records(data):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "stream.pst"
        path.write_bytes(data)
        return list(obdeect_stream.Reader(str(path)))


def _executable():
    exe = find_program("obdeect.propagation", "obdeect-propagation-run")
    if exe is None:
        pytest.skip("obdeect propagation run not installed")
    return exe


def ground_stream(n, seed, telescopes=(1, 2), chunk_rows=65536, curvature_applied=None):
    """G-frame photons with positions, directions, altitudes and integer weights.

    curvature_applied: run key (0 or 1) if not None, absent otherwise.
    """
    rng = np.random.default_rng(seed)
    writer = StreamWriter()
    run = {
        "run_number": 1,
        "observation_level_m": SITE_M,
        "telescope_x_m": [100.0, -50.0],
        "telescope_y_m": [0.0, 300.0],
        "telescope_z_m": [0.0, 0.0],
        "telescope_r_m": [12.0, 12.0],
        "refraction_applied": 0,
        "wavelength_min_nm": 250.0,
        "wavelength_max_nm": 700.0,
    }
    if curvature_applied is not None:
        run["curvature_applied"] = curvature_applied
    writer.metadata("run", run)
    writer.metadata("event", {"event_number": 5})
    for telescope in telescopes:
        writer.metadata("telescope", {"telescope_id": telescope, "frame": 1})
        for first in range(0, n, chunk_rows):
            rows = min(chunk_rows, n - first)
            columns = {
                "x": rng.uniform(-300, 300, rows),
                "y": rng.uniform(-300, 300, rows),
                "cx": rng.uniform(-0.3, 0.3, rows),
                "cy": rng.uniform(-0.3, 0.3, rows),
                "t": rng.uniform(0, 10, rows),
                "zem": rng.uniform(2000, 8000, rows),
                "w": rng.integers(1, 5, rows).astype(float),
            }
            writer.photons(telescope, 1, columns, first_index=first)
        writer.metadata("telescope_end", {"telescope_id": telescope, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    return writer.getvalue()


@pytest.fixture
def envelope_files(tmp_path):
    qe = tmp_path / "qe.dat"
    qe.write_text("# qe x 0.9\n200 0.9\n800 0.5\n", encoding="utf-8")
    mirror = tmp_path / "mirror.dat"
    mirror.write_text("200 0.8\n800 0.8\n", encoding="utf-8")
    return {"qe": qe, "mirror": mirror}


def _efficiency_status(records):
    return next(r["metadata"] for r in records if r["type"] == "efficiency_status")


def test_export_is_lazy():
    assert callable(obdeect_propagation.propagate)
    assert "propagate" in obdeect_propagation.__all__


def test_efficiency_status_and_counts(tmp_path, envelope_files):
    n = 30000
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(n, seed=1))
    thinned_out = tmp_path / "thinned.pst"
    weighted_out = tmp_path / "weighted.pst"
    options = {
        "transmission": TABLE_PATH,
        "site_altitude_m": SITE_M,
        "seed": 11,
        "envelope": envelope_files,
        "store_p": True,
    }
    thinned = obdeect_propagation.propagate(source, thinned_out, mode="thinned", **options)
    weighted = obdeect_propagation.propagate(source, weighted_out, mode="weighted", **options)

    assert thinned["photons_in"] == 2 * n
    assert thinned["dropped"] == 0
    records = read_records(thinned_out.read_bytes())
    efficiency_status = _efficiency_status(records)
    assert efficiency_status["factors"] == "atmosphere,wavelength,qe,mirror"
    assert efficiency_status["atmosphere.state"] == "applied"
    assert efficiency_status["atmosphere.source"] == str(TABLE_PATH)
    assert efficiency_status["wavelength.state"] == "applied"
    assert efficiency_status["qe.state"] == "envelope"
    assert efficiency_status["mirror.state"] == "envelope"
    assert len(efficiency_status["qe.wavelength_nm"]) == len(efficiency_status["qe.value"]) == 2
    # The efficiency status follows the run record.
    assert [r["type"] for r in records[:3]] == ["stream", "run", "efficiency_status"]

    # Thinned count against the weighted expectation: w_out = w p, so the variance of the
    # number of survivors is sum w p (1 - p) = sum w_out (1 - p).
    weighted_records = read_records(weighted_out.read_bytes())
    w_out = np.concatenate([
        r["columns"]["w"] for r in weighted_records if r["type"] == "photons"
    ]).astype(float)
    p = np.concatenate([
        r["columns"]["p"] for r in weighted_records if r["type"] == "photons"
    ]).astype(float)
    mean = w_out.sum()
    sigma = np.sqrt(np.sum(w_out * (1.0 - p)))
    survivors = sum(len(r["columns"]["w"]) for r in records if r["type"] == "photons")
    assert abs(survivors - mean) < 4.0 * sigma
    assert weighted["weight_out"] == pytest.approx(mean, rel=1e-6)


def test_input_efficiency_status_showing_atmosphere_applied_is_rejected(tmp_path):
    writer = StreamWriter()
    writer.metadata(
        "run",
        {
            "observation_level_m": SITE_M,
            "telescope_x_m": [0.0],
            "telescope_y_m": [0.0],
            "telescope_z_m": [0.0],
            "telescope_r_m": [1.0],
            "refraction_applied": 0,
            "wavelength_min_nm": 300.0,
            "wavelength_max_nm": 600.0,
        },
    )
    writer.metadata("efficiency_status", {"factors": "atmosphere", "atmosphere.state": "applied"})
    writer.metadata("event", {"event_number": 1})
    writer.metadata("event_end")
    writer.metadata("run_end")
    source = tmp_path / "tampered.pst"
    source.write_bytes(writer.getvalue())
    with pytest.raises(ValueError, match="already applied"):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
        )


def test_executable_matches_api(tmp_path, envelope_files):
    exe = _executable()
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(5000, seed=2, chunk_rows=2048))
    api_out = tmp_path / "api.pst"
    obdeect_propagation.propagate(
        source,
        api_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=7,
        envelope=envelope_files,
    )
    cli_out = tmp_path / "cli.pst"
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "7",
        "--envelope",
        f"qe={envelope_files['qe']}",
        "--envelope",
        f"mirror={envelope_files['mirror']}",
        "-i",
        str(source),
        "-o",
        str(cli_out),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "photons in 10000" in result.stderr
    assert cli_out.read_bytes() == api_out.read_bytes()


def test_executable_stdin_matches_file_input(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(5000, seed=2, chunk_rows=2048))
    command = [
        _executable(),
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "7",
        "-i",
    ]
    file_result = subprocess.run(
        [*command, str(source), "-o", "-"], capture_output=True, check=False
    )
    stdin_result = subprocess.run(
        [*command, "-", "-o", "-"], input=source.read_bytes(), capture_output=True, check=False
    )
    assert file_result.returncode == 0, file_result.stderr
    assert stdin_result.returncode == 0, stdin_result.stderr
    assert b"photons in 10000" in stdin_result.stderr
    assert stdin_result.stderr == file_result.stderr
    assert stdin_result.stdout == file_result.stdout


SCATTERING_TABLE = DATA / "scattering_exponential.dat"
SCATTER_TELESCOPES = [(1, 3000.0, 0.0, 0.0, 12.0), (2, 0.0, -2500.0, 0.0, 12.0)]


def beam_stream(n, chunk_rows=65536, radii=True):
    """Vertical beam at x = 0, 15 km above the site (sky frame, 355 nm), telescope 1 only.

    Telescope 1 is 3 km from the beam and telescope 2 (2.5 km south) has no input block.
    """
    run = {
        "run_number": 1,
        "observation_level_m": SITE_M,
        "telescope_x_m": [3000.0, 0.0],
        "telescope_y_m": [0.0, -2500.0],
        "telescope_z_m": [0.0, 0.0],
        "refraction_applied": 0,
        "wavelength_min_nm": 250.0,
        "wavelength_max_nm": 700.0,
    }
    if radii:
        run["telescope_r_m"] = [12.0, 12.0]
    writer = StreamWriter()
    writer.metadata("run", run)
    writer.metadata("event", {"event_number": 7})
    writer.metadata("telescope", {"telescope_id": 1, "frame": 0})
    for first in range(0, n, chunk_rows):
        rows = min(chunk_rows, n - first)
        i = np.arange(first, first + rows)
        columns = {
            "x": np.zeros(rows),
            "y": np.zeros(rows),
            "z": np.full(rows, 15000.0),
            "cx": np.zeros(rows),
            "cy": np.zeros(rows),
            "cz": -np.ones(rows),
            "t": np.zeros(rows),
            "w": 1.0 + i % 3,
            "wavelength": np.full(rows, 355.0),
        }
        writer.photons(1, 0, columns, first_index=first)
    writer.metadata("telescope_end", {"telescope_id": 1, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    return writer.getvalue()


@pytest.fixture
def mirror_file(tmp_path):
    mirror = tmp_path / "mirror.dat"
    mirror.write_text("200 0.8\n800 0.8\n", encoding="utf-8")
    return mirror


def test_scattering_adds_kernel_weight_times_envelope(tmp_path, mirror_file):
    n = 400
    source = tmp_path / "beam.pst"
    source.write_bytes(beam_stream(n))
    options = {
        "transmission": TABLE_PATH,
        "site_altitude_m": SITE_M,
        "seed": 1,
        "envelope": {"mirror": mirror_file},
        "mode": "weighted",
    }
    base = obdeect_propagation.propagate(source, tmp_path / "off.pst", **options)
    on_out = tmp_path / "on.pst"
    on = obdeect_propagation.propagate(
        source, on_out, scattering=SCATTERING_TABLE, scattering_points=8, **options
    )

    # Primaries unchanged; the scattered weight is the kernel result (same seed and event, the
    # primary index of the beam photons) times the envelope 0.8. It is not extinguished again.
    i = np.arange(n)
    columns = {
        "x": np.zeros(n),
        "y": np.zeros(n),
        "z": np.full(n, 15000.0),
        "cx": np.zeros(n),
        "cy": np.zeros(n),
        "t": np.zeros(n),
        "w": 1.0 + i % 3,
        "wavelength": np.full(n, 355.0),
    }
    kernel = obdeect_propagation.scatter(
        columns,
        obdeect_propagation.ScatteringTable(str(SCATTERING_TABLE)),
        obdeect_propagation.TransmissionTable(
            str(TABLE_PATH), site_altitude_m=SITE_M, wavelength_interpolation="floor"
        ),
        telescopes=SCATTER_TELESCOPES,
        seed=1,
        event=7,
        first_index=0,
        n_points=8,
    )
    expected = 0.8 * sum(
        float(np.sum(photons["w"][photons["cz"] < 0], dtype=np.float64))
        for photons in kernel.values()
    )
    assert on["weight_out"] - base["weight_out"] == pytest.approx(expected, rel=1e-5)
    assert on["photons_out"] > base["photons_out"]

    records = read_records(on_out.read_bytes())
    assert any(r["type"] == "photons" and r["telescope_id"] == 2 for r in records)
    efficiency_status = _efficiency_status(records)
    assert efficiency_status["factors"] == "atmosphere,wavelength,scattering,mirror"
    assert efficiency_status["scattering.state"] == "applied"
    assert efficiency_status["scattering.source"] == str(SCATTERING_TABLE)
    assert efficiency_status["scattering.points_per_photon"] == 8


def test_scattering_chunk_sizes_agree(tmp_path, mirror_file):
    results = []
    for chunk_rows in (3000, 333):
        source = tmp_path / f"beam_{chunk_rows}.pst"
        source.write_bytes(beam_stream(3000, chunk_rows=chunk_rows))
        results.append(
            obdeect_propagation.propagate(
                source,
                tmp_path / f"out_{chunk_rows}.pst",
                transmission=TABLE_PATH,
                site_altitude_m=SITE_M,
                seed=3,
                envelope={"mirror": mirror_file},
                scattering=SCATTERING_TABLE,
            )
        )
    assert results[0]["photons_out"] == results[1]["photons_out"] > 0
    assert results[0]["weight_out"] == pytest.approx(results[1]["weight_out"], rel=1e-12)


def test_scattering_rejects_invalid_settings(tmp_path):
    source = tmp_path / "beam.pst"
    source.write_bytes(beam_stream(10, radii=False))
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 1}
    with pytest.raises(ValueError, match="telescope_r_m"):
        obdeect_propagation.propagate(
            source, tmp_path / "out.pst", scattering=SCATTERING_TABLE, **options
        )
    source.write_bytes(beam_stream(10))
    with pytest.raises(ValueError):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            scattering=SCATTERING_TABLE,
            scattering_points=0,
            **options,
        )


def test_executable_matches_api_with_scattering(tmp_path):
    exe = _executable()
    source = tmp_path / "beam.pst"
    source.write_bytes(beam_stream(300))
    api_out = tmp_path / "api.pst"
    obdeect_propagation.propagate(
        source,
        api_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=5,
        scattering=SCATTERING_TABLE,
        scattering_points=4,
    )
    cli_out = tmp_path / "cli.pst"
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "5",
        "--scattering",
        str(SCATTERING_TABLE),
        "--scattering-points",
        "4",
        "-i",
        str(source),
        "-o",
        str(cli_out),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert cli_out.read_bytes() == api_out.read_bytes()


def test_rejects_identical_paths(tmp_path):
    source = tmp_path / "same.pst"
    original = ground_stream(100, seed=3, chunk_rows=64)
    source.write_bytes(original)
    alias = tmp_path / ".." / tmp_path.name / "same.pst"
    for output in (source, alias):
        with pytest.raises(ValueError, match="same file"):
            obdeect_propagation.propagate(
                source,
                output,
                transmission=TABLE_PATH,
                site_altitude_m=SITE_M,
                seed=1,
            )
    assert source.read_bytes() == original


def test_executable_rejects_identical_paths(tmp_path):
    exe = _executable()
    source = tmp_path / "same.pst"
    original = ground_stream(100, seed=3, chunk_rows=64)
    source.write_bytes(original)
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "1",
        "-i",
        str(source),
        "-o",
        str(source),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 2, result.stderr
    assert "same file" in result.stderr
    assert source.read_bytes() == original


def _photon_block(records, telescope, name):
    return np.concatenate([
        r["columns"][name]
        for r in records
        if r["type"] == "photons" and r["telescope_id"] == telescope
    ]).astype(float)


def test_pathlen_default_and_override(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(4000, seed=3))
    out = tmp_path / "out.pst"
    obdeect_propagation.propagate(
        source,
        out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=4,
        mode="weighted",
        store_p=True,
        pathlen_m=28.0,
        wavelength_interpolation="floor",
        telescope_pathlen_m={2: 0.0},
    )
    records = read_records(out.read_bytes())
    assert list(_efficiency_status(records)["atmosphere.pathlen_m"]) == [28.0, 0.0]

    # Acceptance of each row with od_focus = tau(h + l) - tau(h).
    table = obdeect_propagation.TransmissionTable(
        str(TABLE_PATH), SITE_M, wavelength_interpolation="floor"
    )
    for telescope, pathlen in ((1, 28.0), (2, 0.0)):
        lam = _photon_block(records, telescope, "wavelength")
        zem = _photon_block(records, telescope, "zem")
        cx = _photon_block(records, telescope, "cx")
        cy = _photon_block(records, telescope, "cy")
        p = _photon_block(records, telescope, "p")
        assert p.size > 0
        od_tel = table.optical_depth(lam, SITE_M)
        od_focus = table.optical_depth(lam, SITE_M + pathlen) - od_tel
        airmass = 1.0 / np.sqrt(1.0 - cx**2 - cy**2)
        expected = table.transmission2(lam, zem, airmass, od_tel, od_focus)
        assert np.abs(p - expected).max() <= 1e-6


def test_pathlen_changes_acceptance(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(4000, seed=3))
    with_path = tmp_path / "with.pst"
    without = tmp_path / "without.pst"
    base = {
        "transmission": TABLE_PATH,
        "site_altitude_m": SITE_M,
        "seed": 4,
        "mode": "weighted",
        "store_p": True,
    }
    obdeect_propagation.propagate(source, with_path, pathlen_m=28.0, **base)
    obdeect_propagation.propagate(source, without, **base)
    p_with = _photon_block(read_records(with_path.read_bytes()), 1, "p")
    p_without = _photon_block(read_records(without.read_bytes()), 1, "p")
    assert p_with.shape == p_without.shape
    assert np.abs(p_with - p_without).max() > 1e-5
    status = _efficiency_status(read_records(without.read_bytes()))
    assert list(status["atmosphere.pathlen_m"]) == [0.0, 0.0]


@pytest.mark.parametrize(
    "options",
    [
        {"pathlen_m": -1.0},
        {"pathlen_m": float("nan")},
        {"pathlen_m": float("inf")},
        {"telescope_pathlen_m": {1: -0.5}},
        {"telescope_pathlen_m": {1: float("nan")}},
        {"telescope_pathlen_m": {0: 1.0}},
        {"telescope_pathlen_m": {3: 1.0}},
    ],
)
def test_pathlen_rejects_invalid(tmp_path, options):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(100, seed=5))
    with pytest.raises(ValueError):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            **options,
        )


def test_pathlen_ids_must_be_integers(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(100, seed=5))
    with pytest.raises(TypeError):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            telescope_pathlen_m={1.5: 1.0},
        )


def test_executable_pathlen_matches_api(tmp_path):
    exe = _executable()
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(5000, seed=6, chunk_rows=2048))
    api_out = tmp_path / "api.pst"
    obdeect_propagation.propagate(
        source,
        api_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=7,
        pathlen_m=28.0,
        telescope_pathlen_m={2: 0.0},
    )
    cli_out = tmp_path / "cli.pst"
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "7",
        "--pathlen",
        "28",
        "--telescope-pathlen",
        "2=0",
        "-i",
        str(source),
        "-o",
        str(cli_out),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert cli_out.read_bytes() == api_out.read_bytes()
    status = _efficiency_status(read_records(cli_out.read_bytes()))
    assert list(status["atmosphere.pathlen_m"]) == [28.0, 0.0]


@pytest.mark.parametrize(
    "args",
    [
        ["--pathlen", "-1"],
        ["--pathlen", "nan"],
        ["--pathlen", "inf"],
        ["--telescope-pathlen", "1"],
        ["--telescope-pathlen", "0=1"],
        ["--telescope-pathlen", "x=1"],
        ["--telescope-pathlen", "1=-1"],
        ["--telescope-pathlen", "1=nan"],
        ["--telescope-pathlen", "1=2", "--telescope-pathlen", "1=3"],
    ],
)
def test_executable_rejects_invalid_pathlen(args):
    exe = _executable()
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "1",
        *args,
    ]
    result = subprocess.run(command, capture_output=True, text=True, input="", check=False)
    assert result.returncode == 2, result.stderr


PROFILE_PATH = DATA / "atmprof_ecmwf_north_winter_fixed.dat"


def sky_stream(n, seed, curvature_applied=None, refraction_applied=0):
    """Sky-frame photons with random positions and downward directions, telescope 1 and 2.

    curvature_applied: run key (0 or 1) if not None, absent otherwise.
    refraction_applied: run key (0: photons not refracted, 1: already refracted).
    """
    rng = np.random.default_rng(seed)
    writer = StreamWriter()
    run = {
        "run_number": 1,
        "observation_level_m": SITE_M,
        "telescope_x_m": [100.0, -50.0],
        "telescope_y_m": [0.0, 300.0],
        "telescope_z_m": [0.0, 0.0],
        "refraction_applied": refraction_applied,
        "wavelength_min_nm": 250.0,
        "wavelength_max_nm": 700.0,
    }
    if curvature_applied is not None:
        run["curvature_applied"] = curvature_applied
    writer.metadata("run", run)
    writer.metadata("event", {"event_number": 5})
    for telescope in (1, 2):
        writer.metadata("telescope", {"telescope_id": telescope, "frame": 0})
        cx = rng.uniform(-0.3, 0.3, n)
        cy = rng.uniform(-0.3, 0.3, n)
        columns = {
            "x": rng.uniform(-3000, 3000, n),
            "y": rng.uniform(-3000, 3000, n),
            "z": rng.uniform(3000, 15000, n),
            "cx": cx,
            "cy": cy,
            "cz": -np.sqrt(1.0 - cx**2 - cy**2),
            "t": rng.uniform(0, 10, n),
            "w": rng.integers(1, 5, n).astype(float),
        }
        writer.photons(telescope, 0, columns)
        writer.metadata("telescope_end", {"telescope_id": telescope, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    return writer.getvalue()


# Float32 directions (cx, cy, cz) of photons with |cz| between 6e-4 and 2e-3. Each is a unit
# vector to float32 precision; an absolute check on cz rejected them. The last one has
# cx^2 + cy^2 rounded above 1.
NEAR_HORIZONTAL_DIRECTIONS = (
    (0.960585952, 0.277977228, -0.001825427),
    (0.946776211, -0.321888804, -0.001533536),
    (0.999714673, 0.023859283, -0.001124879),
    (0.989117742, -0.147124648, -0.000614837),
    (0.975054383, -0.221962571, -0.001266957),
    (1.0000001192, 0.0, -0.0001),
)


def near_horizontal_stream():
    """Sky-frame photons on telescope 1 with the directions of NEAR_HORIZONTAL_DIRECTIONS."""
    cx, cy, cz = (np.array(column, dtype=np.float32) for column in zip(*NEAR_HORIZONTAL_DIRECTIONS))
    n = len(cx)
    rng = np.random.default_rng(11)
    writer = StreamWriter()
    writer.metadata(
        "run",
        {
            "run_number": 1,
            "observation_level_m": SITE_M,
            "telescope_x_m": [0.0],
            "telescope_y_m": [0.0],
            "telescope_z_m": [0.0],
            "refraction_applied": 0,
            "wavelength_min_nm": 250.0,
            "wavelength_max_nm": 700.0,
        },
    )
    writer.metadata("event", {"event_number": 1})
    writer.metadata("telescope", {"telescope_id": 1, "frame": 0})
    columns = {
        "x": rng.uniform(-3000, 3000, n),
        "y": rng.uniform(-3000, 3000, n),
        "z": rng.uniform(3000, 15000, n),
        "cx": cx,
        "cy": cy,
        "cz": cz,
        "t": rng.uniform(0, 10, n),
        "w": np.ones(n),
    }
    writer.photons(1, 0, columns)
    writer.metadata("telescope_end", {"telescope_id": 1, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    return writer.getvalue()


def test_near_horizontal_photons_pass_checks_and_reach_ground(tmp_path):
    # Geometry only: an atmospheric profile or transmission absorbs photons at this grazing angle.
    source = tmp_path / "in.pst"
    source.write_bytes(near_horizontal_stream())
    summary = obdeect_propagation.propagate(
        source, tmp_path / "out.pst", extinction=False, site_altitude_m=SITE_M, seed=1
    )
    assert summary["photons_in"] == len(NEAR_HORIZONTAL_DIRECTIONS)
    assert summary["photons_out"] == len(NEAR_HORIZONTAL_DIRECTIONS)


def _geometry_options(**extra):
    # 400 to 700 nm: the extinction is far from zero, so the normal run drops no row.
    return {
        "site_altitude_m": SITE_M,
        "seed": 3,
        "profile": PROFILE_PATH,
        "wavelength_range": (400.0, 700.0),
        **extra,
    }


def test_geometry_only_matches_geometry_of_normal_run(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(sky_stream(3000, seed=8))
    normal_out = tmp_path / "normal.pst"
    geometry_out = tmp_path / "geometry.pst"
    obdeect_propagation.propagate(
        source,
        normal_out,
        transmission=TABLE_PATH,
        mode="weighted",
        store_p=True,
        **_geometry_options(),
    )
    summary = obdeect_propagation.propagate(
        source, geometry_out, extinction=False, **_geometry_options()
    )
    normal = read_records(normal_out.read_bytes())
    geometry = read_records(geometry_out.read_bytes())

    for telescope in (1, 2):
        assert _photon_block(normal, telescope, "p").min() > 0.0
        for name in ("x", "y", "cx", "cy", "t", "wavelength"):
            np.testing.assert_array_equal(
                _photon_block(geometry, telescope, name), _photon_block(normal, telescope, name)
            )
        # Weights are the input weights (integers 1 to 4), not multiplied by p.
        w = _photon_block(geometry, telescope, "w")
        assert np.all(w == np.round(w)) and w.min() >= 1.0
        assert w.sum() > _photon_block(normal, telescope, "w").sum()
    assert "p" not in next(r for r in geometry if r["type"] == "photons")["columns"]
    assert summary["photons_out"] == summary["photons_in"] - summary["dropped"]
    assert summary["weight_out"] == pytest.approx(
        sum(_photon_block(geometry, t, "w").sum() for t in (1, 2))
    )


def test_geometry_only_efficiency_status_is_pending(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(sky_stream(500, seed=9))
    out = tmp_path / "out.pst"
    obdeect_propagation.propagate(source, out, extinction=False, **_geometry_options())
    efficiency_status = _efficiency_status(read_records(out.read_bytes()))
    assert efficiency_status["factors"] == "atmosphere,wavelength"
    assert efficiency_status["atmosphere.state"] == "pending"
    assert "atmosphere.source" not in efficiency_status
    assert efficiency_status["wavelength.state"] == "applied"


def test_geometry_only_on_ground_stream_keeps_every_photon(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(1000, seed=2))
    out = tmp_path / "out.pst"
    summary = obdeect_propagation.propagate(
        source, out, extinction=False, site_altitude_m=SITE_M, seed=1, mode="weighted"
    )
    assert summary["photons_out"] == summary["photons_in"] == 2000
    records = read_records(out.read_bytes())
    assert _photon_block(records, 1, "wavelength").min() >= 250.0
    original = read_records(source.read_bytes())
    for name in ("x", "y", "cx", "cy", "t", "zem", "w"):
        np.testing.assert_array_equal(
            _photon_block(records, 2, name), _photon_block(original, 2, name)
        )


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"transmission": TABLE_PATH}, "transmission"),
        ({"envelope": {"qe": "qe.dat"}}, "envelope"),
        ({"mode": "thinned"}, "thinned"),
        ({"scattering": "table.dat"}, "scattering"),
        ({"scattering_points": 4}, "scattering_points"),
        # The default value is rejected too, as the CLI rejects any --scattering-points.
        ({"scattering_points": 8}, "scattering_points"),
        ({"pathlen_m": 0.0}, "pathlen_m"),
        ({"store_p": True}, "store_p"),
        ({"pathlen_m": 1.0}, "pathlen_m"),
        ({"telescope_pathlen_m": {1: 1.0}}, "telescope_pathlen_m"),
    ],
)
def test_geometry_only_rejects_extinction_options(tmp_path, options, message):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(10, seed=1))
    with pytest.raises(ValueError, match=message):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            extinction=False,
            site_altitude_m=SITE_M,
            seed=1,
            **options,
        )


def test_transmission_is_required_with_extinction(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(10, seed=1))
    with pytest.raises(ValueError, match="transmission is required"):
        obdeect_propagation.propagate(source, tmp_path / "out.pst", site_altitude_m=SITE_M, seed=1)


def test_executable_geometry_only_matches_api(tmp_path):
    exe = _executable()
    source = tmp_path / "in.pst"
    source.write_bytes(sky_stream(2000, seed=4))
    api_out = tmp_path / "api.pst"
    obdeect_propagation.propagate(source, api_out, extinction=False, **_geometry_options())
    cli_out = tmp_path / "cli.pst"
    command = [
        exe,
        "--no-extinction",
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "3",
        "--profile",
        str(PROFILE_PATH),
        "--wavelength-range",
        "400",
        "700",
        "-i",
        str(source),
        "-o",
        str(cli_out),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert cli_out.read_bytes() == api_out.read_bytes()
    assert _efficiency_status(read_records(cli_out.read_bytes()))["atmosphere.state"] == "pending"


@pytest.mark.parametrize(
    "args",
    [
        ["--transmission", str(TABLE_PATH)],
        ["--envelope", "qe=qe.dat"],
        ["--mode", "thinned"],
        ["--scattering", "table.dat"],
        ["--scattering-points", "8"],
        ["--store-p"],
        ["--pathlen", "1"],
        ["--telescope-pathlen", "1=1"],
    ],
)
def test_executable_geometry_only_rejects_extinction_options(args):
    exe = _executable()
    command = [exe, "--no-extinction", "--site-altitude", str(SITE_M), "--seed", "1", *args]
    result = subprocess.run(command, capture_output=True, text=True, input="", check=False)
    assert result.returncode == 2, result.stderr
    assert "--no-extinction" in result.stderr


def test_executable_geometry_only_accepts_weighted_mode(tmp_path):
    exe = _executable()
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(100, seed=1))
    command = [
        exe,
        "--no-extinction",
        "--mode",
        "weighted",
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "1",
        "-i",
        str(source),
        "-o",
        str(tmp_path / "out.pst"),
    ]
    assert subprocess.run(command, capture_output=True, check=False).returncode == 0


def test_curvature_default_refuses_corsika_curved_input(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(100, seed=1, curvature_applied=1))
    out = tmp_path / "out.pst"
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 1}
    with pytest.raises(ValueError, match="CORSIKA CURVED"):
        obdeect_propagation.propagate(source, out, **options)
    summary = obdeect_propagation.propagate(source, out, curvature="corsika", **options)
    assert summary["photons_in"] == 200


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"curvature": "curved"}, "curvature must be one of"),
        ({"curvature": "Corsika"}, "curvature must be one of"),
        ({"curvature": "corsika", "earth_radius_m": 1e5}, "above 1e6 m"),
        ({"curvature": "spherical", "earth_radius_m": float("nan")}, "above 1e6 m"),
    ],
)
def test_curvature_rejects_invalid_settings(tmp_path, kwargs, message):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(10, seed=1, curvature_applied=1))
    with pytest.raises(ValueError, match=message):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            **kwargs,
        )


def test_spherical_needs_uncurved_input(tmp_path):
    source = tmp_path / "in.pst"
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 1}
    source.write_bytes(sky_stream(10, seed=1, curvature_applied=1))
    with pytest.raises(ValueError, match="already applied"):
        obdeect_propagation.propagate(
            source, tmp_path / "out.pst", curvature="spherical", **options
        )
    source.write_bytes(sky_stream(10, seed=1))
    with pytest.raises(ValueError, match="curvature_applied"):
        obdeect_propagation.propagate(
            source, tmp_path / "out.pst", curvature="spherical", **options
        )


def test_spherical_accepts_profile(tmp_path):
    # Spherical refraction and air delay: the run keeps refraction_applied = 0 (bending).
    source = tmp_path / "in.pst"
    source.write_bytes(sky_stream(10, seed=1, curvature_applied=0))
    out = tmp_path / "out.pst"
    summary = obdeect_propagation.propagate(
        source,
        out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=1,
        curvature="spherical",
        profile=PROFILE_PATH,
    )
    assert summary["photons_in"] == 20
    records = read_records(out.read_bytes())
    run = next(r["metadata"] for r in records if r["type"] == "run")
    assert run["curvature_applied"] == 1
    assert run["refraction_applied"] == 1
    assert _efficiency_status(records)["atmosphere.geometry"] == "spherical"


def test_spherical_status_and_run_record(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(sky_stream(500, seed=3, curvature_applied=0))
    out = tmp_path / "out.pst"
    obdeect_propagation.propagate(
        source,
        out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=3,
        curvature="spherical",
        earth_radius_m=6.0e6,
    )
    records = read_records(out.read_bytes())
    status = _efficiency_status(records)
    assert status["atmosphere.geometry"] == "spherical"
    assert status["atmosphere.earth_radius_m"] == pytest.approx(6.0e6)
    run = next(r["metadata"] for r in records if r["type"] == "run")
    assert run["curvature_applied"] == 1


def test_corsika_status_and_no_status_without_curvature(tmp_path):
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 2}
    corsika_in = tmp_path / "corsika.pst"
    corsika_in.write_bytes(ground_stream(100, seed=2, curvature_applied=1))
    corsika_out = tmp_path / "corsika_out.pst"
    obdeect_propagation.propagate(corsika_in, corsika_out, curvature="corsika", **options)
    status = _efficiency_status(read_records(corsika_out.read_bytes()))
    assert status["atmosphere.geometry"] == "corsika-curved"
    assert status["atmosphere.earth_radius_m"] == pytest.approx(6371315.0)

    plain_in = tmp_path / "plain.pst"
    plain_in.write_bytes(ground_stream(100, seed=2))
    plain_out = tmp_path / "plain_out.pst"
    obdeect_propagation.propagate(plain_in, plain_out, **options)
    status = _efficiency_status(read_records(plain_out.read_bytes()))
    assert "atmosphere.geometry" not in status
    assert "atmosphere.earth_radius_m" not in status


def test_spherical_geometry_only_is_allowed(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(sky_stream(200, seed=4, curvature_applied=0))
    summary = obdeect_propagation.propagate(
        source,
        tmp_path / "out.pst",
        site_altitude_m=SITE_M,
        seed=4,
        extinction=False,
        curvature="spherical",
    )
    assert summary["photons_in"] == 400


def test_executable_curvature_matches_api(tmp_path):
    exe = _executable()
    source = tmp_path / "in.pst"
    source.write_bytes(sky_stream(2000, seed=5, curvature_applied=0))
    api_out = tmp_path / "api.pst"
    obdeect_propagation.propagate(
        source,
        api_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=7,
        curvature="spherical",
        earth_radius_m=6.0e6,
    )
    cli_out = tmp_path / "cli.pst"
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "7",
        "--curvature",
        "spherical",
        "--earth-radius-m",
        "6000000",
        "-i",
        str(source),
        "-o",
        str(cli_out),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert cli_out.read_bytes() == api_out.read_bytes()


def test_executable_rejects_unknown_curvature(tmp_path):
    exe = _executable()
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "1",
        "--curvature",
        "flat",
    ]
    result = subprocess.run(command, capture_output=True, text=True, input="", check=False)
    assert result.returncode == 2, result.stderr
    assert "--curvature" in result.stderr


def _without_telescopes(records, dropped):
    """The records of a run without the blocks of the telescopes in `dropped`.

    A block is its telescope record, photons, telescope status and telescope_end; run-level
    records are kept.
    """
    kept, block = [], None
    for record in records:
        kind = record["type"]
        if kind == "telescope":
            block = record["metadata"]["telescope_id"]
        if kind == "photons":
            drop = record["telescope_id"] in dropped
        elif block is not None and kind in ("telescope", "efficiency_status", "telescope_end"):
            drop = block in dropped
        else:
            drop = False
        if not drop:
            kept.append(record)
        if kind == "telescope_end":
            block = None
    return kept


def test_telescopes_keeps_selected_blocks_and_run_records(tmp_path):
    n = 3000
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(n, seed=4, chunk_rows=1024))
    full_out = tmp_path / "full.pst"
    selected_out = tmp_path / "selected.pst"
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 9}
    obdeect_propagation.propagate(source, full_out, **options)
    summary = obdeect_propagation.propagate(source, selected_out, telescopes=[2], **options)

    full = read_records(full_out.read_bytes())
    selected = read_records(selected_out.read_bytes())
    assert summary["photons_in"] == n  # the rows of telescope 2 only
    expected = _without_telescopes(full, dropped={1})
    assert [r["type"] for r in selected] == [r["type"] for r in expected]
    assert {r["metadata"]["telescope_id"] for r in selected if r["type"] == "telescope"} == {2}
    for got, want in zip(selected, expected, strict=True):
        if got["type"] == "photons":
            assert got["telescope_id"] == want["telescope_id"]
            assert got["first_index"] == want["first_index"]
            for name in want["columns"]:
                np.testing.assert_array_equal(got["columns"][name], want["columns"][name])
        else:
            assert repr(got.get("metadata")) == repr(want.get("metadata"))


def test_telescopes_executable_matches_api(tmp_path):
    exe = _executable()
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(2000, seed=6, chunk_rows=512))
    api_out = tmp_path / "api.pst"
    obdeect_propagation.propagate(
        source,
        api_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=3,
        telescopes=[2, 1],
    )
    cli_out = tmp_path / "cli.pst"
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "3",
        "--telescopes",
        "2,1,2",
        "-i",
        str(source),
        "-o",
        str(cli_out),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "photons in 4000" in result.stderr
    assert cli_out.read_bytes() == api_out.read_bytes()


@pytest.mark.parametrize("telescopes", [[], [0], [3], [1.5]])
def test_telescopes_rejects_invalid_selection(tmp_path, telescopes):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(100, seed=8))
    with pytest.raises((ValueError, TypeError)):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            telescopes=telescopes,
        )


def test_telescopes_combine_with_scattering(tmp_path):
    """A telescope selection with scattering keeps the scattering target and drops the others.

    The beam's primaries are on telescope 1, the scattered photons reach telescope 2. Keeping
    telescope 2 gives its block of the full run, unchanged.
    """
    source = tmp_path / "beam.pst"
    source.write_bytes(beam_stream(400))
    options = {
        "transmission": TABLE_PATH,
        "site_altitude_m": SITE_M,
        "seed": 1,
        "mode": "weighted",
        "wavelength_interpolation": "floor",
    }
    full_out = tmp_path / "full.pst"
    obdeect_propagation.propagate(source, full_out, scattering=SCATTERING_TABLE, **options)
    kept_out = tmp_path / "kept.pst"
    kept = obdeect_propagation.propagate(
        source, kept_out, scattering=SCATTERING_TABLE, telescopes=[2], **options
    )

    full = read_records(full_out.read_bytes())
    kept_records = read_records(kept_out.read_bytes())
    telescopes = {r["metadata"]["telescope_id"] for r in kept_records if r["type"] == "telescope"}
    assert telescopes == {2}
    full_weight_2 = sum(
        float(np.sum(r["columns"]["w"], dtype=np.float64))
        for r in full
        if r["type"] == "photons" and r["telescope_id"] == 2
    )
    assert full_weight_2 > 0.0
    assert kept["weight_out"] == pytest.approx(full_weight_2, rel=1e-6)


@pytest.mark.parametrize("selection", ["0", "1,,2", "x", "2,"])
def test_executable_rejects_invalid_telescopes(selection):
    exe = _executable()
    command = [
        exe,
        "--transmission",
        str(TABLE_PATH),
        "--site-altitude",
        str(SITE_M),
        "--seed",
        "1",
        "--telescopes",
        selection,
    ]
    result = subprocess.run(command, capture_output=True, text=True, input="", check=False)
    assert result.returncode == 2, result.stderr
    assert "--telescopes" in result.stderr


def test_threads_output_is_identical_to_single_thread(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(3000, seed=4, chunk_rows=1024))
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 9}
    single = tmp_path / "single.pst"
    multi = tmp_path / "multi.pst"
    obdeect_propagation.propagate(source, single, **options)
    obdeect_propagation.propagate(source, multi, threads=4, **options)
    assert multi.read_bytes() == single.read_bytes()


@pytest.mark.parametrize("threads", [0, -1, 257, 1.5])
def test_threads_rejects_invalid_count(tmp_path, threads):
    source = tmp_path / "in.pst"
    source.write_bytes(ground_stream(100, seed=8))
    with pytest.raises((ValueError, TypeError)):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            threads=threads,
        )


def _ground_wavelength_stream(wavelengths, direction_cx=0.2, zem_m=2300.0):
    """G-frame photons of telescope 1, one per wavelength (nm), weight 1."""
    n = len(wavelengths)
    writer = StreamWriter()
    writer.metadata("run", {**_geometry_run(), "refraction_applied": 0})
    writer.metadata("event", {"event_number": 1})
    writer.metadata("telescope", {"telescope_id": 1, "frame": 1})
    writer.photons(
        1,
        1,
        {
            "x": np.zeros(n),
            "y": np.zeros(n),
            "cx": np.full(n, direction_cx),
            "cy": np.zeros(n),
            "t": np.zeros(n),
            "zem": np.full(n, zem_m),
            "w": np.ones(n),
            "wavelength": np.asarray(wavelengths, dtype=float),
        },
    )
    writer.metadata("telescope_end", {"telescope_id": 1, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    return writer.getvalue()


def _geometry_run():
    return {
        "run_number": 1,
        "observation_level_m": SITE_M,
        "telescope_x_m": [0.0],
        "telescope_y_m": [0.0],
        "telescope_z_m": [0.0],
        "wavelength_min_nm": 250.0,
        "wavelength_max_nm": 700.0,
    }


def test_wavelength_interpolation_default_is_linear_and_is_recorded(tmp_path):
    source = tmp_path / "in.pst"
    source.write_bytes(_ground_wavelength_stream([400.25, 400.999, 401.0]))
    default_out = tmp_path / "default.pst"
    obdeect_propagation.propagate(
        source,
        default_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=1,
        mode="weighted",
    )
    floor_out = tmp_path / "floor.pst"
    obdeect_propagation.propagate(
        source,
        floor_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=1,
        mode="weighted",
        wavelength_interpolation="floor",
    )
    default_status = _efficiency_status(read_records(default_out.read_bytes()))
    floor_status = _efficiency_status(read_records(floor_out.read_bytes()))
    assert default_status["atmosphere.wavelength_interpolation"] == "linear"
    assert floor_status["atmosphere.wavelength_interpolation"] == "floor"


def test_linear_mode_is_continuous_and_interpolates_the_envelope(tmp_path):
    """Linear mode: w(400.999) is close to w(401.0); floor mode uses the row of 400 nm for both.

    The envelope factor is linear in wavelength between its nodes (390 nm 0.9, 410 nm 0.5), so
    E(400.25) = 0.9 - 0.4 * 10.25 / 20 in linear mode, E(400) = 0.7 in floor mode.
    """
    wavelengths = [400.0, 400.25, 400.999, 401.0]
    source = tmp_path / "in.pst"
    source.write_bytes(_ground_wavelength_stream(wavelengths))
    envelope = tmp_path / "envelope.dat"
    envelope.write_text("250 0.9\n700 0.5\n", encoding="utf-8")
    results = {}
    for mode in ("linear", "floor"):
        out = tmp_path / f"{mode}.pst"
        obdeect_propagation.propagate(
            source,
            out,
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            mode="weighted",
            envelope={"eff": envelope},
            pathlen_m=28.0,
            wavelength_interpolation=mode,
        )
        records = read_records(out.read_bytes())
        results[mode] = _photon_block(records, 1, "w")
    linear, floor = results["linear"], results["floor"]
    assert floor[1] == floor[0] and floor[2] == floor[0]
    assert abs(linear[2] - linear[3]) < 1e-3
    assert linear[1] != linear[0]
    assert linear[0] == pytest.approx(floor[0], rel=1e-6)  # integer nm: the modes agree

    # Against the table with the envelope at the linear wavelength (the stage's arithmetic).
    table = obdeect_propagation.TransmissionTable(
        str(TABLE_PATH), SITE_M, wavelength_interpolation="linear"
    )
    depths = table.telescope_depths(SITE_M, 28.0)
    lam = np.asarray(wavelengths)
    od_tel = np.interp(lam, np.arange(len(depths[0])), depths[0])
    od_focus = np.interp(lam, np.arange(len(depths[1])), depths[1])
    airmass = np.full(lam.shape, 1.0 / np.sqrt(1.0 - 0.2**2))
    e = np.interp(lam, [250.0, 700.0], [0.9, 0.5])
    expected = table.transmission2(lam, np.full(lam.shape, 2300.0), airmass, od_tel, od_focus) * e
    np.testing.assert_allclose(linear, expected, rtol=1e-5)


def test_dispersion_needs_profile_and_unrefracted_input(tmp_path):
    source = tmp_path / "sky.pst"
    source.write_bytes(sky_stream(50, seed=2))
    with pytest.raises(ValueError, match="profile"):
        obdeect_propagation.propagate(
            source,
            tmp_path / "a.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            dispersion_reference_nm=500.0,
        )
    refracted = tmp_path / "refracted.pst"
    refracted.write_bytes(sky_stream(50, seed=2, refraction_applied=1))
    with pytest.raises(ValueError, match="refract"):
        obdeect_propagation.propagate(
            refracted,
            tmp_path / "b.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            profile=PROFILE_PATH,
            dispersion_reference_nm=500.0,
        )
    with pytest.raises(ValueError):
        obdeect_propagation.propagate(
            source,
            tmp_path / "c.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            profile=PROFILE_PATH,
            dispersion_reference_nm=150.0,
        )


def test_dispersion_changes_bent_photons_by_a_small_amount(tmp_path):
    """Dispersion scales the bending by n(lambda)-1 relative to the 500 nm profile.

    Photons at 250 to 700 nm. The position changes by the phase factor k - 1 of the refraction
    (about 8 % of it at 250 nm); the largest change seen is about 2 cm. The time changes by the
    group factor k_group - 1 of the air delay, which is larger in the UV (up to about 1.6 ns for
    the photons of this stream, whose air delay is several ns). The status records the value.
    """
    source = tmp_path / "sky.pst"
    source.write_bytes(sky_stream(200, seed=3))
    options = dict(
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=4,
        mode="weighted",
        profile=PROFILE_PATH,
    )
    off_out = tmp_path / "off.pst"
    obdeect_propagation.propagate(source, off_out, **options)
    on_out = tmp_path / "on.pst"
    obdeect_propagation.propagate(source, on_out, dispersion_reference_nm=500.0, **options)
    off = read_records(off_out.read_bytes())
    on = read_records(on_out.read_bytes())
    assert "atmosphere.dispersion_reference_nm" not in _efficiency_status(off)
    assert _efficiency_status(on)["atmosphere.dispersion_reference_nm"] == 500.0
    for column, bound in (("x", 0.05), ("y", 0.05), ("t", 2.0)):  # m for x and y, ns for t
        delta = _photon_block(on, 1, column) - _photon_block(off, 1, column)
        assert np.any(delta != 0.0)
        assert np.max(np.abs(delta)) < bound


def test_near_horizon_photons_are_bent_not_dropped(tmp_path):
    """Zenith 88 degrees, beyond the series limit of the bending: bent by quadrature."""
    writer = StreamWriter()
    writer.metadata(
        "run",
        {**_geometry_run(), "refraction_applied": 0},
    )
    writer.metadata("event", {"event_number": 1})
    writer.metadata("telescope", {"telescope_id": 1, "frame": 0})
    cx = np.array([np.sin(np.radians(88.0))])
    cx = cx.astype(np.float32).astype(float)
    cz = -np.sqrt(1.0 - cx**2)
    writer.photons(
        1,
        0,
        {
            "x": np.zeros(1),
            "y": np.zeros(1),
            "z": np.full(1, 10000.0),
            "cx": cx,
            "cy": np.zeros(1),
            "cz": cz,
            "t": np.zeros(1),
            "w": np.ones(1),
        },
    )
    writer.metadata("telescope_end", {"telescope_id": 1, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    source = tmp_path / "horizon.pst"
    source.write_bytes(writer.getvalue())
    summary = obdeect_propagation.propagate(
        source,
        tmp_path / "out.pst",
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=1,
        mode="weighted",
        profile=PROFILE_PATH,
    )
    assert summary["dropped"] == 0
    assert summary["photons_out"] == 1


def _refracted_ground_stream(n, seed):
    """G-frame photons of two telescopes, refraction already applied (run key 1)."""
    rng = np.random.default_rng(seed)
    writer = StreamWriter()
    run = _geometry_run()
    run.update({
        "telescope_x_m": [100.0, -50.0],
        "telescope_y_m": [0.0, 300.0],
        "telescope_z_m": [0.0, 0.0],
        "refraction_applied": 1,
    })
    writer.metadata("run", run)
    writer.metadata("event", {"event_number": 5})
    for telescope in (1, 2):
        writer.metadata("telescope", {"telescope_id": telescope, "frame": 1})
        cx = rng.uniform(-0.25, 0.25, n)
        cy = rng.uniform(-0.25, 0.25, n)
        writer.photons(
            telescope,
            1,
            {
                "x": rng.uniform(-3000, 3000, n),
                "y": rng.uniform(-3000, 3000, n),
                "cx": cx,
                "cy": cy,
                "t": rng.uniform(0, 1000, n),
                "zem": rng.uniform(2200, 12000, n),
                "w": rng.integers(1, 5, n).astype(float),
            },
        )
        writer.metadata("telescope_end", {"telescope_id": telescope, "photons": 0.0})
    writer.metadata("event_end")
    writer.metadata("run_end")
    return writer.getvalue()


def test_to_sky_round_trip_through_the_stage(tmp_path):
    """to_sky, then the stage without extinction, gives back the G-frame photons."""
    n = 300
    ground = tmp_path / "ground.pst"
    ground.write_bytes(_refracted_ground_stream(n, seed=6))
    sky = tmp_path / "sky.pst"
    converted = obdeect_propagation.to_sky(ground, sky, profile=PROFILE_PATH)
    assert converted == {"photons_in": 2 * n, "photons_out": 2 * n, "dropped": 0}
    back = tmp_path / "back.pst"
    summary = obdeect_propagation.propagate(
        sky, back, extinction=False, site_altitude_m=SITE_M, seed=1, profile=PROFILE_PATH
    )
    assert summary["photons_out"] == 2 * n
    original = read_records(ground.read_bytes())
    returned = read_records(back.read_bytes())
    for telescope in (1, 2):
        for column in ("x", "y", "cx", "cy", "t", "zem"):
            want = _photon_block(original, telescope, column)
            got = _photon_block(returned, telescope, column)
            np.testing.assert_allclose(got, want, rtol=1e-6, atol=1e-2, err_msg=column)


def test_to_sky_rejects_sky_frame_input(tmp_path):
    source = tmp_path / "sky.pst"
    source.write_bytes(sky_stream(10, seed=1))
    with pytest.raises(ValueError, match="frame G"):
        obdeect_propagation.to_sky(source, tmp_path / "out.pst")
    with pytest.raises(ValueError):
        obdeect_propagation.to_sky(source, source)


def test_geometry_only_rejects_wavelength_interpolation(tmp_path):
    source = tmp_path / "sky.pst"
    source.write_bytes(sky_stream(10, seed=1))
    with pytest.raises(ValueError, match="wavelength_interpolation"):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out.pst",
            extinction=False,
            site_altitude_m=SITE_M,
            seed=1,
            wavelength_interpolation="floor",
        )
    with pytest.raises(ValueError, match="wavelength_interpolation must be one of"):
        obdeect_propagation.propagate(
            source,
            tmp_path / "out2.pst",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            wavelength_interpolation="nearest",
        )
