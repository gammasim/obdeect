"""Photon-stream round trips through the Python API and the obdeect stream cat / tee tools."""

import subprocess

import numpy as np
import obdeect.stream as obdeect_stream
import pytest
from obdeect.stream.cli import find_program

REQUIRED = ("x", "y", "cx", "cy", "t", "w")


def _columns(rows, names, offset=0):
    return {
        name: np.arange(rows, dtype=np.float32) * np.float32(0.5) + np.float32(i + offset)
        for i, name in enumerate(names)
    }


def _write_stream(path, rows=1000, telescope_id=1):
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {"run_number": 1})
        writer.write("telescope", {"telescope_id": telescope_id, "frame": "telescope"})
        writer.write_photons(telescope_id, "telescope", 0, _columns(rows, REQUIRED))
        writer.write("telescope_end", {"photons": float(rows)})
        writer.write("run_end", {})


def test_round_trip(tmp_path):
    path = tmp_path / "a.pst"
    first = _columns(1000, REQUIRED + ("zem", "wavelength"))
    second = _columns(7, REQUIRED, offset=100)
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write(
            "run",
            {
                "run_number": 3,
                "telescope_x_m": np.array([1.0, 2.0, -3.5]),
                "name": "site",
                "energy_min_gev": 0.25,
                "refraction_applied": True,
            },
        )
        writer.write("telescope", {"telescope_id": 1, "frame": 2})
        writer.write_photons(1, "telescope", 0, first)
        writer.write_photons(1, 2, 1000, second)
        writer.write_blob(1204, 7, b"\x01\x02\x03")
        writer.write("telescope_end", {"photons": 1007.0})
        writer.write("run_end", {})

    records = list(obdeect_stream.Reader(str(path)))
    assert [r["type"] for r in records] == [
        "stream",
        "run",
        "telescope",
        "photons",
        "photons",
        "blob",
        "telescope_end",
        "run_end",
    ]
    assert records[0]["metadata"]["producer"] == "pytest"
    run = records[1]["metadata"]
    assert run["run_number"] == 3
    assert run["name"] == "site"
    assert run["energy_min_gev"] == 0.25
    assert run["refraction_applied"] == 1
    np.testing.assert_array_equal(run["telescope_x_m"], [1.0, 2.0, -3.5])
    assert run["telescope_x_m"].dtype == np.float64

    for record, expected in ((records[3], first), (records[4], second)):
        assert record["telescope_id"] == 1
        assert set(record["columns"]) == set(expected)
        assert record["first_index"] in (0, 1000)
        for name, values in expected.items():
            got = record["columns"][name]
            assert got.dtype == np.float32
            assert np.array_equal(got, values)
    assert records[4]["frame"] == 2
    assert records[4]["first_index"] == 1000

    blob = records[5]
    assert blob["eventio_type"] == 1204
    assert blob["eventio_id"] == 7
    assert blob["data"] == b"\x01\x02\x03"


def test_empty_chunk(tmp_path):
    path = tmp_path / "empty.pst"
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write_photons(5, "camera", 42, {})
        writer.write_photons(5, "camera", 42, {name: np.empty(0, np.float32) for name in REQUIRED})
    chunks = [r for r in obdeect_stream.Reader(str(path)) if r["type"] == "photons"]
    assert len(chunks) == 2
    for chunk in chunks:
        assert chunk["columns"] == {}
        assert chunk["frame"] == 3
        assert chunk["first_index"] == 42


def test_invalid_input_rejected(tmp_path):
    with obdeect_stream.Writer(str(tmp_path / "bad.pst"), "pytest") as writer:
        with pytest.raises(ValueError):
            writer.write_photons(1, "ground", 0, {"x": np.zeros(2, np.float32)})
        with pytest.raises(ValueError):
            writer.write_photons(
                1, "ground", 0, _columns(2, REQUIRED) | {"x": np.zeros(3, np.float32)}
            )
        with pytest.raises(TypeError):
            writer.write_photons(1, "ground", 0, _columns(2, REQUIRED) | {"x": np.zeros(2)})
        with pytest.raises(ValueError):
            writer.write_photons(
                1, "ground", 0, _columns(2, REQUIRED) | {"q": np.zeros(2, np.float32)}
            )
        with pytest.raises(ValueError):
            writer.write("photons", {})
        with pytest.raises(ValueError):
            writer.write_photons(1, "moon", 0, {})


def test_garbage_input_rejected(tmp_path):
    path = tmp_path / "garbage.pst"
    path.write_bytes(b"this is not a photon stream at all")
    with pytest.raises(RuntimeError):
        list(obdeect_stream.Reader(str(path)))


def test_missing_stream_end_rejected(tmp_path):
    path = tmp_path / "full.pst"
    _write_stream(path)
    data = path.read_bytes()
    cut = tmp_path / "cut.pst"
    cut.write_bytes(data[:-16])  # drop the 16-byte stream_end record
    with pytest.raises(RuntimeError):
        list(obdeect_stream.Reader(str(cut)))
    empty = tmp_path / "empty.pst"
    empty.write_bytes(b"")
    with pytest.raises(RuntimeError):
        list(obdeect_stream.Reader(str(empty)))


def test_data_after_stream_end_rejected(tmp_path):
    path = tmp_path / "full.pst"
    _write_stream(path)
    data = path.read_bytes()
    for name, extra in (("one_byte.pst", b"\x00"), ("concatenated.pst", data)):
        bad = tmp_path / name
        bad.write_bytes(data + extra)
        with pytest.raises(RuntimeError, match="after stream_end"):
            list(obdeect_stream.Reader(str(bad)))


def test_non_integral_frame_rejected(tmp_path):
    with obdeect_stream.Writer(str(tmp_path / "frame.pst"), "pytest") as writer:
        for frame in (2.9, True):
            with pytest.raises(TypeError):
                writer.write_photons(1, frame, 0, _columns(2, REQUIRED))


def test_write_after_close_rejected(tmp_path):
    writer = obdeect_stream.Writer(str(tmp_path / "closed.pst"), "pytest")
    writer.close()
    writer.close()
    with pytest.raises(RuntimeError):
        writer.write("run", {})


def test_abort_is_terminal(tmp_path):
    writer = obdeect_stream.Writer(str(tmp_path / "aborted.pst"), "pytest")
    writer.abort()
    for operation in (lambda: writer.write("run", {}), writer.close):
        with pytest.raises(RuntimeError, match="closed"):
            operation()


@pytest.mark.skipif(
    find_program("obdeect.stream", "obdeect-stream-cat") is None,
    reason="obdeect stream cat not installed",
)
def test_ps_cat(tmp_path):
    cat = find_program("obdeect.stream", "obdeect-stream-cat")
    a, b = tmp_path / "a.pst", tmp_path / "b.pst"
    _write_stream(a, rows=1000, telescope_id=1)
    _write_stream(b, rows=500, telescope_id=2)

    proc = subprocess.run([cat, str(a), str(b)], capture_output=True, check=True)
    out = tmp_path / "joined.pst"
    out.write_bytes(proc.stdout)
    records = list(obdeect_stream.Reader(str(out)))
    types = [r["type"] for r in records]
    assert types.count("stream") == 1
    assert types.count("run") == 2
    photons = {r["telescope_id"]: len(r["columns"]["w"]) for r in records if r["type"] == "photons"}
    assert photons == {1: 1000, 2: 500}

    stdin = subprocess.run([cat, "-"], input=a.read_bytes(), capture_output=True, check=True)
    assert [r["type"] for r in obdeect_stream.Reader(_bytes_path(tmp_path, stdin.stdout))].count(
        "photons"
    ) == 1

    summary = subprocess.run([cat, "--summary", str(a), str(b)], capture_output=True, check=True)
    text = summary.stdout.decode()
    assert "total: 1500 photons" in text
    assert "telescope 1: 1000 photons" in text
    assert "telescope 2: 500 photons" in text

    truncated = tmp_path / "truncated.pst"
    truncated.write_bytes(a.read_bytes()[:-20])
    failed = subprocess.run([cat, str(truncated)], capture_output=True)
    assert failed.returncode != 0


@pytest.mark.skipif(
    find_program("obdeect.stream", "obdeect-stream-tee") is None,
    reason="obdeect stream tee not installed",
)
def test_ps_tee(tmp_path):
    tee = find_program("obdeect.stream", "obdeect-stream-tee")
    source = tmp_path / "source.pst"
    _write_stream(source)
    copy = tmp_path / "copy.pst"
    proc = subprocess.run(
        [tee, str(copy)], input=source.read_bytes(), capture_output=True, check=True
    )
    assert proc.stdout == source.read_bytes()
    assert copy.read_bytes() == source.read_bytes()


@pytest.mark.skipif(
    find_program("obdeect.stream", "obdeect-stream-tee") is None,
    reason="obdeect stream tee not installed",
)
@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("missing_header", "does not start with a stream record"),
        ("missing_stream_end", "no stream_end"),
        ("trailing_data", "after stream_end"),
        ("empty", "empty input"),
    ],
)
def test_ps_tee_rejects_bad_record_order(tmp_path, case, message):
    tee = find_program("obdeect.stream", "obdeect-stream-tee")
    source = tmp_path / "source.pst"
    _write_stream(source)
    data = source.read_bytes()
    header_size = 16 + int.from_bytes(data[8:16], "little")  # stream record: header + payload
    bad = {
        "missing_header": data[header_size:],
        "missing_stream_end": data[:-16],
        "trailing_data": data + data[-16:],
        "empty": b"",
    }[case]
    copy = tmp_path / "copy.pst"
    proc = subprocess.run([tee, str(copy)], input=bad, capture_output=True)
    assert proc.returncode != 0
    assert message in proc.stderr.decode()


def _bytes_path(directory, data):
    path = directory / "stdout.pst"
    path.write_bytes(data)
    return str(path)
