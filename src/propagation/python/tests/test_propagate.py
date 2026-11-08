"""Atmosphere stage through the Python API and the dlufd-propagate executable."""

import io
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import dlufd
import numpy as np
import pytest

try:
    import gstrickd
except ImportError:
    gstrickd = None

# gstrickd's own Reader and Writer are used when importable; the local codec below is only a
# fallback. TODO: remove the local codec once gstrickd >= 0.1 is on PyPI.
USE_GSTRICKD_CODEC = all(hasattr(gstrickd, name) for name in ("Reader", "Writer"))

DATA = Path(__file__).resolve().parents[2] / "cpp" / "tests" / "data"
TABLE_PATH = DATA / "atm_trans_2000_1_10_0_0_2000.dat"
SITE_M = 2000.0
MAGIC = 0x52545347
RECORD_TYPES = {
    "stream": 1,
    "run": 2,
    "event": 3,
    "telescope": 4,
    "photons": 5,
    "telescope_end": 6,
    "event_end": 7,
    "run_end": 8,
    "blob": 9,
    "ledger": 10,
    "stream_end": 15,
}
NAMES = {value: key for key, value in RECORD_TYPES.items()}
COLUMNS = ["x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "wavelength", "p"]


# Local writer and reader of the gstrickd stream format (docs/format.md of gstrickd); fallback only.
def _metadata(values):
    out = [struct.pack("<I", len(values))]
    for key, value in values.items():
        encoded = key.encode()
        out.append(struct.pack("<H", len(encoded)) + encoded)
        if isinstance(value, int):
            out.append(struct.pack("<Bq", 1, value))
        elif isinstance(value, float):
            out.append(struct.pack("<Bd", 2, value))
        elif isinstance(value, str):
            raw = value.encode()
            out.append(struct.pack("<BI", 3, len(raw)) + raw)
        else:
            array = np.asarray(value, dtype="<f8")
            out.append(struct.pack("<BI", 4, array.size) + array.tobytes())
    return b"".join(out)


def _record(kind, payload):
    return struct.pack("<IHHQ", MAGIC, RECORD_TYPES[kind], 0, len(payload)) + payload


class _LocalStreamWriter:
    def __init__(self):
        self.buffer = io.BytesIO()
        self.buffer.write(
            _record(
                "stream",
                _metadata({"format_major": 1, "format_minor": 0, "producer": "test"}),
            )
        )

    def metadata(self, kind, values=None):
        self.buffer.write(_record(kind, _metadata(values or {})))

    def photons(self, telescope_id, frame, columns, first_index=0):
        rows = len(next(iter(columns.values())))
        mask = sum(1 << COLUMNS.index(name) for name in columns)
        header = struct.pack("<IIIIQ", telescope_id, frame, mask, rows, first_index)
        body = b"".join(
            np.asarray(columns[name], dtype="<f4").tobytes() for name in COLUMNS if name in columns
        )
        self.buffer.write(_record("photons", header + body))

    def getvalue(self):
        self.buffer.write(_record("stream_end", b""))
        return self.buffer.getvalue()


def _read_records_local(data):
    records = []
    pos = 0
    while pos < len(data):
        magic, kind, _, length = struct.unpack_from("<IHHQ", data, pos)
        assert magic == MAGIC
        payload = data[pos + 16 : pos + 16 + length]
        pos += 16 + length
        name = NAMES[kind]
        if name == "stream_end":
            break
        if name == "photons":
            telescope_id, frame, mask, rows, first_index = struct.unpack_from("<IIIIQ", payload)
            offset = 24
            columns = {}
            for bit, column in enumerate(COLUMNS):
                if mask >> bit & 1:
                    columns[column] = np.frombuffer(payload, "<f4", rows, offset)
                    offset += 4 * rows
            records.append({
                "type": name,
                "telescope_id": telescope_id,
                "frame": frame,
                "first_index": first_index,
                "columns": columns,
            })
        elif name in ("blob", "stream"):
            records.append({"type": name})
        else:
            records.append({"type": name, "metadata": _decode_metadata(payload)})
    return records


def _decode_metadata(payload):
    (count,) = struct.unpack_from("<I", payload, 0)
    pos = 4
    values = {}
    for _ in range(count):
        (key_size,) = struct.unpack_from("<H", payload, pos)
        key = payload[pos + 2 : pos + 2 + key_size].decode()
        pos += 2 + key_size
        (kind,) = struct.unpack_from("<B", payload, pos)
        pos += 1
        if kind == 1:
            (value,) = struct.unpack_from("<q", payload, pos)
            pos += 8
        elif kind == 2:
            (value,) = struct.unpack_from("<d", payload, pos)
            pos += 8
        elif kind == 3:
            (size,) = struct.unpack_from("<I", payload, pos)
            value = payload[pos + 4 : pos + 4 + size].decode()
            pos += 4 + size
        else:
            (size,) = struct.unpack_from("<I", payload, pos)
            value = np.frombuffer(payload, "<f8", size, pos + 4).tolist()
            pos += 4 + 8 * size
        values[key] = value
    return values


class _GstrickdStreamWriter:
    """Writes through gstrickd.Writer into a temporary file; getvalue() returns its bytes."""

    def __init__(self):
        self._directory = tempfile.TemporaryDirectory()
        self._path = Path(self._directory.name) / "stream.ps"
        self._writer = gstrickd.Writer(str(self._path), "test")

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


def _read_records_gstrickd(data):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "stream.ps"
        path.write_bytes(data)
        return list(gstrickd.Reader(str(path)))


StreamWriter = _GstrickdStreamWriter if USE_GSTRICKD_CODEC else _LocalStreamWriter


def read_records(data):
    return _read_records_gstrickd(data) if USE_GSTRICKD_CODEC else _read_records_local(data)


def ground_stream(n, seed, telescopes=(1, 2), chunk_rows=65536):
    """G-frame photons with positions, directions, altitudes and integer weights."""
    rng = np.random.default_rng(seed)
    writer = StreamWriter()
    writer.metadata(
        "run",
        {
            "run_number": 1,
            "observation_level_m": SITE_M,
            "telescope_x_m": [100.0, -50.0],
            "telescope_y_m": [0.0, 300.0],
            "telescope_z_m": [0.0, 0.0],
            "telescope_r_m": [12.0, 12.0],
            "refraction_applied": 0,
            "wavelength_min_nm": 250.0,
            "wavelength_max_nm": 700.0,
        },
    )
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


def _ledger(records):
    return next(r["metadata"] for r in records if r["type"] == "ledger")


def test_export_is_lazy():
    assert callable(dlufd.propagate)
    assert "propagate" in dlufd.__all__


def test_ledger_and_counts(tmp_path, envelope_files):
    n = 30000
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(n, seed=1))
    thinned_out = tmp_path / "thinned.ps"
    weighted_out = tmp_path / "weighted.ps"
    options = {
        "transmission": TABLE_PATH,
        "site_altitude_m": SITE_M,
        "seed": 11,
        "envelope": envelope_files,
        "store_p": True,
    }
    thinned = dlufd.propagate(source, thinned_out, mode="thinned", **options)
    weighted = dlufd.propagate(source, weighted_out, mode="weighted", **options)

    assert thinned["photons_in"] == 2 * n
    assert thinned["dropped"] == 0
    records = read_records(thinned_out.read_bytes())
    ledger = _ledger(records)
    assert ledger["factors"] == "atmosphere,wavelength,qe,mirror"
    assert ledger["atmosphere.state"] == "applied"
    assert ledger["atmosphere.source"] == str(TABLE_PATH)
    assert ledger["wavelength.state"] == "applied"
    assert ledger["qe.state"] == "envelope"
    assert ledger["mirror.state"] == "envelope"
    assert len(ledger["qe.wavelength_nm"]) == len(ledger["qe.value"]) == 2
    # The ledger follows the run record.
    assert [r["type"] for r in records[:3]] == ["stream", "run", "ledger"]

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


def test_input_ledger_showing_atmosphere_applied_is_rejected(tmp_path):
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
    writer.metadata("ledger", {"factors": "atmosphere", "atmosphere.state": "applied"})
    writer.metadata("event", {"event_number": 1})
    writer.metadata("event_end")
    writer.metadata("run_end")
    source = tmp_path / "tampered.ps"
    source.write_bytes(writer.getvalue())
    with pytest.raises(ValueError, match="already applied"):
        dlufd.propagate(
            source,
            tmp_path / "out.ps",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
        )


def test_executable_matches_api(tmp_path, envelope_files):
    exe = shutil.which("dlufd-propagate") or str(Path(sys.executable).with_name("dlufd-propagate"))
    if not Path(exe).exists():
        pytest.skip("dlufd-propagate not installed")
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(5000, seed=2, chunk_rows=2048))
    api_out = tmp_path / "api.ps"
    dlufd.propagate(
        source,
        api_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=7,
        envelope=envelope_files,
    )
    cli_out = tmp_path / "cli.ps"
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
    source = tmp_path / "beam.ps"
    source.write_bytes(beam_stream(n))
    options = {
        "transmission": TABLE_PATH,
        "site_altitude_m": SITE_M,
        "seed": 1,
        "envelope": {"mirror": mirror_file},
        "mode": "weighted",
    }
    base = dlufd.propagate(source, tmp_path / "off.ps", **options)
    on_out = tmp_path / "on.ps"
    on = dlufd.propagate(
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
    kernel = dlufd.scatter(
        columns,
        dlufd.ScatteringTable(str(SCATTERING_TABLE)),
        dlufd.TransmissionTable(str(TABLE_PATH), site_altitude_m=SITE_M),
        telescopes=SCATTER_TELESCOPES,
        seed=1,
        event=7,
        first_index=0,
        n_points=8,
    )
    expected = 0.8 * sum(
        float(np.sum(photons["w"], dtype=np.float64)) for photons in kernel.values()
    )
    assert on["weight_out"] - base["weight_out"] == pytest.approx(expected, rel=1e-5)
    assert on["photons_out"] > base["photons_out"]

    records = read_records(on_out.read_bytes())
    assert any(r["type"] == "photons" and r["telescope_id"] == 2 for r in records)
    ledger = _ledger(records)
    assert ledger["factors"] == "atmosphere,wavelength,scattering,mirror"
    assert ledger["scattering.state"] == "applied"
    assert ledger["scattering.source"] == str(SCATTERING_TABLE)
    assert ledger["scattering.points_per_photon"] == 8


def test_scattering_chunk_sizes_agree(tmp_path, mirror_file):
    results = []
    for chunk_rows in (3000, 333):
        source = tmp_path / f"beam_{chunk_rows}.ps"
        source.write_bytes(beam_stream(3000, chunk_rows=chunk_rows))
        results.append(
            dlufd.propagate(
                source,
                tmp_path / f"out_{chunk_rows}.ps",
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
    source = tmp_path / "beam.ps"
    source.write_bytes(beam_stream(10, radii=False))
    options = {"transmission": TABLE_PATH, "site_altitude_m": SITE_M, "seed": 1}
    with pytest.raises(ValueError, match="telescope_r_m"):
        dlufd.propagate(source, tmp_path / "out.ps", scattering=SCATTERING_TABLE, **options)
    source.write_bytes(beam_stream(10))
    with pytest.raises(ValueError):
        dlufd.propagate(
            source, tmp_path / "out.ps", scattering=SCATTERING_TABLE, scattering_points=0, **options
        )


def test_executable_matches_api_with_scattering(tmp_path):
    exe = shutil.which("dlufd-propagate") or str(Path(sys.executable).with_name("dlufd-propagate"))
    if not Path(exe).exists():
        pytest.skip("dlufd-propagate not installed")
    source = tmp_path / "beam.ps"
    source.write_bytes(beam_stream(300))
    api_out = tmp_path / "api.ps"
    dlufd.propagate(
        source,
        api_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=5,
        scattering=SCATTERING_TABLE,
        scattering_points=4,
    )
    cli_out = tmp_path / "cli.ps"
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
    source = tmp_path / "same.ps"
    original = ground_stream(100, seed=3, chunk_rows=64)
    source.write_bytes(original)
    alias = tmp_path / ".." / tmp_path.name / "same.ps"
    for output in (source, alias):
        with pytest.raises(ValueError, match="same file"):
            dlufd.propagate(
                source,
                output,
                transmission=TABLE_PATH,
                site_altitude_m=SITE_M,
                seed=1,
            )
    assert source.read_bytes() == original


def test_executable_rejects_identical_paths(tmp_path):
    exe = shutil.which("dlufd-propagate") or str(Path(sys.executable).with_name("dlufd-propagate"))
    if not Path(exe).exists():
        pytest.skip("dlufd-propagate not installed")
    source = tmp_path / "same.ps"
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


def _executable():
    exe = shutil.which("dlufd-propagate") or str(Path(sys.executable).with_name("dlufd-propagate"))
    if not Path(exe).exists():
        pytest.skip("dlufd-propagate not installed")
    return exe


def _photon_block(records, telescope, name):
    return np.concatenate([
        r["columns"][name]
        for r in records
        if r["type"] == "photons" and r["telescope_id"] == telescope
    ]).astype(float)


def test_pathlen_default_and_override(tmp_path):
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(4000, seed=3))
    out = tmp_path / "out.ps"
    dlufd.propagate(
        source,
        out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=4,
        mode="weighted",
        store_p=True,
        pathlen_m=28.0,
        telescope_pathlen_m={2: 0.0},
    )
    records = read_records(out.read_bytes())
    assert list(_ledger(records)["atmosphere.pathlen_m"]) == [28.0, 0.0]

    # Acceptance of each row with od_focus = tau(h + l) - tau(h).
    table = dlufd.TransmissionTable(str(TABLE_PATH), SITE_M)
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
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(4000, seed=3))
    with_path = tmp_path / "with.ps"
    without = tmp_path / "without.ps"
    base = {
        "transmission": TABLE_PATH,
        "site_altitude_m": SITE_M,
        "seed": 4,
        "mode": "weighted",
        "store_p": True,
    }
    dlufd.propagate(source, with_path, pathlen_m=28.0, **base)
    dlufd.propagate(source, without, **base)
    p_with = _photon_block(read_records(with_path.read_bytes()), 1, "p")
    p_without = _photon_block(read_records(without.read_bytes()), 1, "p")
    assert p_with.shape == p_without.shape
    assert np.abs(p_with - p_without).max() > 1e-5
    assert list(_ledger(read_records(without.read_bytes()))["atmosphere.pathlen_m"]) == [0.0, 0.0]


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
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(100, seed=5))
    with pytest.raises(ValueError):
        dlufd.propagate(
            source,
            tmp_path / "out.ps",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            **options,
        )


def test_pathlen_ids_must_be_integers(tmp_path):
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(100, seed=5))
    with pytest.raises(TypeError):
        dlufd.propagate(
            source,
            tmp_path / "out.ps",
            transmission=TABLE_PATH,
            site_altitude_m=SITE_M,
            seed=1,
            telescope_pathlen_m={1.5: 1.0},
        )


def test_executable_pathlen_matches_api(tmp_path):
    exe = _executable()
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(5000, seed=6, chunk_rows=2048))
    api_out = tmp_path / "api.ps"
    dlufd.propagate(
        source,
        api_out,
        transmission=TABLE_PATH,
        site_altitude_m=SITE_M,
        seed=7,
        pathlen_m=28.0,
        telescope_pathlen_m={2: 0.0},
    )
    cli_out = tmp_path / "cli.ps"
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
    assert list(_ledger(read_records(cli_out.read_bytes()))["atmosphere.pathlen_m"]) == [28.0, 0.0]


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


def sky_stream(n, seed):
    """Sky-frame photons with random positions and downward directions, telescope 1 and 2."""
    rng = np.random.default_rng(seed)
    writer = StreamWriter()
    writer.metadata(
        "run",
        {
            "run_number": 1,
            "observation_level_m": SITE_M,
            "telescope_x_m": [100.0, -50.0],
            "telescope_y_m": [0.0, 300.0],
            "telescope_z_m": [0.0, 0.0],
            "refraction_applied": 0,
            "wavelength_min_nm": 250.0,
            "wavelength_max_nm": 700.0,
        },
    )
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
    source = tmp_path / "in.ps"
    source.write_bytes(sky_stream(3000, seed=8))
    normal_out = tmp_path / "normal.ps"
    geometry_out = tmp_path / "geometry.ps"
    dlufd.propagate(
        source,
        normal_out,
        transmission=TABLE_PATH,
        mode="weighted",
        store_p=True,
        **_geometry_options(),
    )
    summary = dlufd.propagate(source, geometry_out, extinction=False, **_geometry_options())
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


def test_geometry_only_ledger_is_pending(tmp_path):
    source = tmp_path / "in.ps"
    source.write_bytes(sky_stream(500, seed=9))
    out = tmp_path / "out.ps"
    dlufd.propagate(source, out, extinction=False, **_geometry_options())
    ledger = _ledger(read_records(out.read_bytes()))
    assert ledger["factors"] == "atmosphere,wavelength"
    assert ledger["atmosphere.state"] == "pending"
    assert "atmosphere.source" not in ledger
    assert ledger["wavelength.state"] == "applied"


def test_geometry_only_on_ground_stream_keeps_every_photon(tmp_path):
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(1000, seed=2))
    out = tmp_path / "out.ps"
    summary = dlufd.propagate(
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
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(10, seed=1))
    with pytest.raises(ValueError, match=message):
        dlufd.propagate(
            source,
            tmp_path / "out.ps",
            extinction=False,
            site_altitude_m=SITE_M,
            seed=1,
            **options,
        )


def test_transmission_is_required_with_extinction(tmp_path):
    source = tmp_path / "in.ps"
    source.write_bytes(ground_stream(10, seed=1))
    with pytest.raises(ValueError, match="transmission is required"):
        dlufd.propagate(source, tmp_path / "out.ps", site_altitude_m=SITE_M, seed=1)


def test_executable_geometry_only_matches_api(tmp_path):
    exe = _executable()
    source = tmp_path / "in.ps"
    source.write_bytes(sky_stream(2000, seed=4))
    api_out = tmp_path / "api.ps"
    dlufd.propagate(source, api_out, extinction=False, **_geometry_options())
    cli_out = tmp_path / "cli.ps"
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
    assert _ledger(read_records(cli_out.read_bytes()))["atmosphere.state"] == "pending"


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
    source = tmp_path / "in.ps"
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
        str(tmp_path / "out.ps"),
    ]
    assert subprocess.run(command, capture_output=True, check=False).returncode == 0
