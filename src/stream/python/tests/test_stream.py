"""Photon-stream round trips through the Python API and the obdeect stream cat / tee tools."""

import struct
import subprocess
from pathlib import Path

import numpy as np
import obdeect.stream as obdeect_stream
import pytest
from obdeect.stream.cli import find_program

REQUIRED = ("x", "y", "cx", "cy", "t", "w")


def _columns(rows, names, offset=0):
    """Photon columns; the directions are unit vectors (cz = -sqrt(1 - cx^2 - cy^2) if present)."""
    index = np.arange(rows, dtype=np.float32)
    columns = {
        name: index * np.float32(0.5) + np.float32(i + offset) for i, name in enumerate(names)
    }
    if "cx" in columns:
        columns["cx"] = index * np.float32(0.0005)
    if "cy" in columns:
        columns["cy"] = index * np.float32(0.0005)
    if "cz" in columns:
        cx = columns["cx"].astype(np.float64)
        cy = columns["cy"].astype(np.float64)
        columns["cz"] = (-np.sqrt(1.0 - cx * cx - cy * cy)).astype(np.float32)
    return columns


def _write_stream(path, rows=1000, telescope_id=1):
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {"run_number": 1})
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": telescope_id, "frame": "telescope"})
        writer.write_photons(telescope_id, "telescope", 0, _columns(rows, REQUIRED))
        writer.write("telescope_end", {"photons": float(rows)})
        writer.write("event_end", {})
        writer.write("run_end", {})


def _write_block(writer, telescope_id, frame, columns, first=0):
    """Run-level records around one telescope block: event, telescope, photons, telescope_end."""
    writer.write("event", {"event_number": 1})
    writer.write("telescope", {"telescope_id": telescope_id, "frame": frame})
    writer.write_photons(telescope_id, frame, first, columns)
    writer.write("telescope_end", {"telescope_id": telescope_id, "photons": 0.0})
    writer.write("event_end", {})


RECORD_TYPE = {
    "stream": 1,
    "run": 2,
    "event": 3,
    "telescope": 4,
    "photons": 5,
    "telescope_end": 6,
    "event_end": 7,
    "run_end": 8,
    "efficiency_status": 10,
    "stream_end": 15,
}


def _encode_metadata(metadata):
    """Metadata payload: u32 count, then per entry u16 key length, key, u8 type and the value."""
    out = struct.pack("<I", len(metadata))
    for key, value in metadata.items():
        name = key.encode()
        out += struct.pack("<H", len(name)) + name
        if isinstance(value, str):
            text = value.encode()
            out += b"\x03" + struct.pack("<I", len(text)) + text
        elif isinstance(value, float):
            out += b"\x02" + struct.pack("<d", value)
        else:
            out += b"\x01" + struct.pack("<q", value)
    return out


def _record(name, metadata=None):
    """One record as raw bytes (magic, type, reserved, length, payload)."""
    payload = b"" if name == "stream_end" else _encode_metadata(metadata or {})
    return struct.pack("<IHHQ", 0x52545347, RECORD_TYPE[name], 0, len(payload)) + payload


def _split_records(data):
    """(type, record bytes) of each record of a stream's bytes."""
    records, at = [], 0
    while at < len(data):
        size = int.from_bytes(data[at + 8 : at + 16], "little")
        records.append((int.from_bytes(data[at + 4 : at + 6], "little"), data[at : at + 16 + size]))
        at += 16 + size
    return records


def _raw_stream(path, records):
    """A stream from raw records (header and stream_end added). The Writer refuses records out of
    order, so the streams that break the nesting are written with these."""
    header = _record(
        "stream",
        {"format_major": obdeect_stream.format_major, "format_minor": 1, "producer": "pytest"},
    )
    Path(path).write_bytes(header + b"".join(records) + _record("stream_end"))
    return path


@pytest.fixture
def photons_record(tmp_path):
    """The bytes of a valid photons record (telescope 1, ground frame, two photons)."""
    path = tmp_path / "photons-source.pst"
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {})
        _write_block(writer, 1, "ground", _columns(2, REQUIRED))
        writer.write("run_end", {})
    (record,) = [
        raw for kind, raw in _split_records(path.read_bytes()) if kind == RECORD_TYPE["photons"]
    ]
    return record


def _set_format_minor(path, minor):
    """Overwrite the format_minor of the stream record (the width of the value does not change)."""
    data = bytearray(Path(path).read_bytes())
    name = b"format_minor"
    at = data.index(name) + len(name)
    assert data[at] == 1  # integer entry
    data[at + 1 : at + 9] = struct.pack("<q", minor)
    Path(path).write_bytes(bytes(data))


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
                "curvature_applied": 1,
            },
        )
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": 2})
        writer.write_photons(1, "telescope", 0, first)
        writer.write_photons(1, 2, 1000, second)
        writer.write_blob(1204, 7, b"\x01\x02\x03")
        writer.write("telescope_end", {"photons": 1007.0})
        writer.write("event_end", {})
        writer.write("run_end", {})

    records = list(obdeect_stream.Reader(str(path)))
    assert [r["type"] for r in records] == [
        "stream",
        "run",
        "event",
        "telescope",
        "photons",
        "photons",
        "blob",
        "telescope_end",
        "event_end",
        "run_end",
    ]
    assert records[0]["metadata"]["producer"] == "pytest"
    run = records[1]["metadata"]
    assert run["run_number"] == 3
    assert run["name"] == "site"
    assert run["energy_min_gev"] == 0.25
    assert run["refraction_applied"] == 1
    assert run["curvature_applied"] == 1
    np.testing.assert_array_equal(run["telescope_x_m"], [1.0, 2.0, -3.5])
    assert run["telescope_x_m"].dtype == np.float64

    for record, expected in ((records[4], first), (records[5], second)):
        assert record["telescope_id"] == 1
        assert set(record["columns"]) == set(expected)
        assert record["first_index"] in (0, 1000)
        for name, values in expected.items():
            got = record["columns"][name]
            assert got.dtype == np.float32
            assert np.array_equal(got, values)
    assert records[5]["frame"] == 2
    assert records[5]["first_index"] == 1000

    blob = records[6]
    assert blob["eventio_type"] == 1204
    assert blob["eventio_id"] == 7
    assert blob["data"] == b"\x01\x02\x03"


def test_empty_chunk(tmp_path):
    path = tmp_path / "empty.pst"
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {})
        writer.write("event", {})
        writer.write("telescope", {"telescope_id": 5, "frame": "camera"})
        writer.write_photons(5, "camera", 42, {})
        writer.write_photons(5, "camera", 42, {name: np.empty(0, np.float32) for name in REQUIRED})
        writer.write("telescope_end", {"telescope_id": 5, "photons": 0.0})
        writer.write("event_end", {})
        writer.write("run_end", {})
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
    find_program("obdeect.stream", "obdeect-stream-cat") is None,
    reason="obdeect stream cat not installed",
)
@pytest.mark.parametrize(
    "case",
    ["photons_outside_telescope", "telescope_end_without_telescope", "run_end_without_run"],
)
def test_ps_cat_rejects_nesting(tmp_path, case, photons_record):
    cat = find_program("obdeect.stream", "obdeect-stream-cat")
    bad = tmp_path / "bad.pst"
    if case == "photons_outside_telescope":
        middle = [_record("run"), _record("event"), photons_record, _record("event_end")]
    elif case == "telescope_end_without_telescope":
        middle = [
            _record("run"),
            _record("event"),
            _record("telescope_end", {"telescope_id": 1, "photons": 0.0}),
            _record("event_end"),
        ]
    else:
        middle = [_record("run"), _record("run_end")]
    _raw_stream(bad, middle + [_record("run_end")])
    failed = subprocess.run([cat, str(bad)], capture_output=True)
    assert failed.returncode != 0
    assert b"malformed stream" in failed.stderr


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


def test_reader_skips_photons_of_other_telescopes(tmp_path):
    path = tmp_path / "two.pst"
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {"run_number": 1})
        writer.write("event", {"event_number": 1})
        for telescope_id in (1, 2):
            writer.write("telescope", {"telescope_id": telescope_id, "frame": "telescope"})
            writer.write_photons(telescope_id, "telescope", 0, _columns(50, REQUIRED))
            writer.write("telescope_end", {"photons": 50.0})
        writer.write("event_end", {})
        writer.write("run_end", {})

    chunks = [r for r in obdeect_stream.Reader(str(path), [2]) if r["type"] == "photons"]
    assert [c["telescope_id"] for c in chunks] == [1, 2]
    assert chunks[0]["columns"] == {}
    assert chunks[0]["first_index"] == 0
    assert set(chunks[1]["columns"]) == set(REQUIRED)
    assert len(chunks[1]["columns"]["w"]) == 50
    everything = [c for c in obdeect_stream.Reader(str(path)) if c["type"] == "photons"]
    np.testing.assert_array_equal(everything[1]["columns"]["x"], chunks[1]["columns"]["x"])


@pytest.mark.parametrize("frame", ["sky", "ground", "telescope", "camera"])
def test_path_column_round_trip(tmp_path, frame):
    path = tmp_path / "path.pst"
    columns = _columns(50, REQUIRED + ("path",))
    columns["path"][0] = np.float32(0.0)
    columns["path"][1] = np.float32(1e-38)
    columns["path"][2] = np.float32(3.4e38)
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {})
        _write_block(writer, 2, frame, columns)
        writer.write("run_end", {})
    (chunk,) = [r for r in obdeect_stream.Reader(str(path)) if r["type"] == "photons"]
    assert set(chunk["columns"]) == set(columns)
    for name, values in columns.items():
        assert np.array_equal(chunk["columns"][name], values)


def test_path_column_must_be_finite_and_not_negative(tmp_path):
    with obdeect_stream.Writer(str(tmp_path / "bad.pst"), "pytest") as writer:
        for bad in (-1.0, np.nan, np.inf):
            columns = _columns(2, REQUIRED + ("path",))
            columns["path"][1] = np.float32(bad)
            with pytest.raises(ValueError):
                writer.write_photons(1, "ground", 0, columns)


def test_ground_representation_default_and_values():
    assert obdeect_stream.ground_representation({}) == "plane"
    assert obdeect_stream.ground_representation({"run_number": 1}) == "plane"
    assert obdeect_stream.ground_representation({"ground_representation": "plane"}) == "plane"
    assert obdeect_stream.ground_representation({"ground_representation": "sphere"}) == "sphere"
    for bad in ("flat", "Sphere", "", 1, None):
        with pytest.raises(ValueError):
            obdeect_stream.ground_representation({"ground_representation": bad})


def test_sphere_run_round_trip(tmp_path):
    path = tmp_path / "sphere.pst"
    rows = _columns(4, ("x", "y", "z", "cx", "cy", "cz", "t", "zem", "w"))
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {"ground_representation": "sphere", "telescope_r_m": np.array([12.0])})
        _write_block(writer, 1, "ground", rows)
        writer.write("run_end", {})
    records = list(obdeect_stream.Reader(str(path)))
    run = next(r for r in records if r["type"] == "run")["metadata"]
    assert obdeect_stream.ground_representation(run) == "sphere"
    (chunk,) = [r for r in records if r["type"] == "photons"]
    assert set(chunk["columns"]) == set(rows)
    assert np.array_equal(chunk["columns"]["z"], rows["z"])


# Reader rejections (docs/stream/format.md): nesting, values, and the efficiency status rules.


@pytest.mark.parametrize(
    "records",
    [
        [("run", {}), ("event", {}), ("photons", 1, "telescope", 0)],
        [("run", {}), ("event", {}), ("telescope_end", {"telescope_id": 1})],
        [("run", {}), ("run_end", {}), ("run_end", {})],
        [("run", {}), ("event", {}), ("event", {})],
        [("event", {})],
        [("run", {}), ("telescope", {"telescope_id": 1})],
        [("run", {}), ("run", {})],
        [("run", {}), ("event", {}), ("run_end", {})],
        [("run", {}), ("event", {}), ("telescope", {"telescope_id": 1}), ("event_end", {})],
    ],
)
def test_reader_rejects_nesting_violations(tmp_path, records, photons_record):
    path = tmp_path / "nesting.pst"
    raw = []
    for record in records:
        if record[0] == "photons":
            raw.append(photons_record)
        else:
            raw.append(_record(*record))
    raw += [_record("event_end"), _record("run_end")]
    _raw_stream(path, raw)
    with pytest.raises(RuntimeError, match="malformed stream"):
        list(obdeect_stream.Reader(str(path)))


def test_reader_accepts_nested_records(tmp_path):
    path = tmp_path / "nested.pst"
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {})
        writer.write_blob(1204, 0, b"\x00")  # a blob may sit between run and run_end
        writer.write("efficiency_status", {"factors": "qe", "qe.state": "pending"})
        _write_block(writer, 1, "ground", _columns(3, REQUIRED))
        writer.write("run_end", {})
    assert [r["type"] for r in obdeect_stream.Reader(str(path))][-1] == "run_end"


@pytest.mark.parametrize(
    "statuses",
    [
        [{"qe": "pending"}, {"qe": "pending"}],
        [{"qe": "applied"}, {"qe": "pending"}],
        [{"qe": "envelope"}, {"qe": "envelope"}],
        [{"qe": "envelope"}, {"qe": "pending"}],
    ],
)
def test_reader_rejects_backward_or_repeated_efficiency_status(tmp_path, statuses):
    path = tmp_path / "eff.pst"
    raw = [_record("run")]
    for status in statuses:
        ((factor, state),) = status.items()
        raw.append(_record("efficiency_status", {"factors": factor, f"{factor}.state": state}))
    raw += [_record("event"), _record("event_end"), _record("run_end")]
    _raw_stream(path, raw)
    with pytest.raises(RuntimeError, match="invalid transition"):
        list(obdeect_stream.Reader(str(path)))


def test_reader_efficiency_status_order_and_scopes(tmp_path):
    path = tmp_path / "eff_order.pst"
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {})
        writer.write("efficiency_status", {"factors": "qe", "qe.state": "pending"})
        writer.write("efficiency_status", {"factors": "qe", "qe.state": "applied"})
        writer.write("event", {})
        writer.write("telescope", {"telescope_id": 1, "frame": "ground"})
        writer.write("efficiency_status", {"factors": "qe", "qe.state": "pending"})
        writer.write_photons(1, "ground", 0, _columns(2, REQUIRED))
        writer.write("telescope_end", {"photons": 2.0})
        writer.write("event_end", {})
        writer.write("run_end", {})
    assert [r["type"] for r in obdeect_stream.Reader(str(path))].count("photons") == 1

    late = _raw_stream(
        tmp_path / "late.pst",
        [
            _record("run"),
            _record("event"),
            _record("efficiency_status", {"factors": "qe", "qe.state": "pending"}),
            _record("event_end"),
            _record("run_end"),
        ],
    )
    with pytest.raises(RuntimeError, match="after the first event"):
        list(obdeect_stream.Reader(str(late)))


@pytest.mark.parametrize(
    ("column", "value"),
    [("x", np.nan), ("y", np.inf), ("t", np.nan), ("w", np.nan), ("zem", np.nan)],
)
def test_reader_rejects_non_finite_values(tmp_path, column, value):
    path = tmp_path / "nan.pst"
    columns = _columns(3, REQUIRED + ("zem",))
    columns[column][1] = np.float32(value)
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {})
        _write_block(writer, 1, "ground", columns)
        writer.write("run_end", {})
    with pytest.raises(RuntimeError, match="not finite"):
        list(obdeect_stream.Reader(str(path)))


def test_reader_rejects_non_unit_direction(tmp_path):
    path = tmp_path / "nonunit.pst"
    columns = _columns(2, REQUIRED + ("cz",))
    columns["cz"][0] = np.float32(-0.5)  # cx, cy are small: cz is not consistent
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {})
        _write_block(writer, 1, "ground", columns)
        writer.write("run_end", {})
    with pytest.raises(RuntimeError, match="not a unit vector"):
        list(obdeect_stream.Reader(str(path)))


def test_reader_upward_photons_by_format_minor(tmp_path):
    columns = _columns(2, REQUIRED + ("cz",))
    columns["cz"] = -columns["cz"]  # upward
    path = tmp_path / "upward.pst"
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {})
        _write_block(writer, 1, "ground", columns)
        writer.write("run_end", {})
    _set_format_minor(path, 0)  # format 1.0 has no upward photons
    with pytest.raises(RuntimeError, match="upward photon"):
        list(obdeect_stream.Reader(str(path)))
    _set_format_minor(path, 1)  # format 1.1 has them
    assert [r["type"] for r in obdeect_stream.Reader(str(path))].count("photons") == 1


def test_reader_failed_state(tmp_path):
    # event inside an event (the Writer refuses it, so the records are raw)
    path = _raw_stream(
        tmp_path / "failed.pst",
        [
            _record("run"),
            _record("event"),
            _record("event"),
            _record("event_end"),
            _record("run_end"),
        ],
    )
    reader = obdeect_stream.Reader(str(path))
    assert next(reader)["type"] == "stream"
    assert next(reader)["type"] == "run"
    assert next(reader)["type"] == "event"
    with pytest.raises(RuntimeError, match="malformed stream"):
        next(reader)
    with pytest.raises(RuntimeError, match="failed earlier"):
        next(reader)


def test_writer_refuses_out_of_order_records(tmp_path):
    path = tmp_path / "refused.pst"
    writer = obdeect_stream.Writer(str(path), "pytest")
    with pytest.raises(ValueError, match="event outside a run"):
        writer.write("event", {})
    writer.write("run", {})
    with pytest.raises(ValueError, match="stream_end inside a run"):
        writer.close()  # the run is still open: refused, and the stream stays open
    writer.write("run_end", {})
    writer.close()
    # The refused records wrote nothing: the stream reads back in order.
    assert [r["type"] for r in obdeect_stream.Reader(str(path))] == ["stream", "run", "run_end"]


def test_reader_accepts_path_like_and_context_manager(tmp_path):
    path = tmp_path / "pathlike.pst"
    _write_stream(path)
    with obdeect_stream.Reader(path) as reader:
        types = [r["type"] for r in reader]
    assert types[0] == "stream" and types[-1] == "run_end"
    with obdeect_stream.Writer(tmp_path / "pathlike-out.pst", "pytest") as writer:
        writer.write("run", {})
        writer.write("run_end", {})
    assert (tmp_path / "pathlike-out.pst").exists()
    reader = obdeect_stream.Reader(str(path))
    reader.close()
    with pytest.raises(ValueError, match="closed"):
        next(reader)


def test_metadata_record_names_come_from_the_core():
    assert obdeect_stream.METADATA_RECORDS == (
        "run",
        "event",
        "telescope",
        "telescope_end",
        "event_end",
        "run_end",
        "efficiency_status",
    )
    from obdeect.stream import stream as stream_module

    assert obdeect_stream.GROUND_REPRESENTATIONS == ("plane", "sphere")
    assert obdeect_stream.GROUND_REPRESENTATIONS is stream_module.GROUND_REPRESENTATIONS
    assert set(stream_module.__all__) <= set(dir(stream_module))
    for name in stream_module.__all__:
        assert getattr(obdeect_stream, name) is getattr(stream_module, name)
