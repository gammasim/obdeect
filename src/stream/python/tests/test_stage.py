"""obdeect.stream.stage.transform and the conformance streams."""

import json
import re
from pathlib import Path

import numpy as np
import obdeect.stream as obdeect_stream
import pytest
from obdeect.stream.stage import efficiency_status_record, transform

CONFORMANCE = Path(__file__).resolve().parents[2] / "cpp" / "tests" / "data" / "conformance"
REQUIRED = ("x", "y", "cx", "cy", "t", "w")


def _chunk(rows, offset=0.0):
    return {
        name: np.arange(rows, dtype=np.float32) + np.float32(i + offset)
        for i, name in enumerate(REQUIRED)
    }


def _write(path, *, efficiency_status=None):
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", {"run_number": 5})
        if efficiency_status:
            writer.write("efficiency_status", efficiency_status)
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 4, "frame": "ground"})
        writer.write_photons(4, "ground", 0, _chunk(6))
        writer.write_photons(4, "ground", 6, _chunk(2))
        writer.write_blob(1204, 1, b"\x01\x02")
        writer.write("telescope_end", {"photons": 8.0})
        writer.write("event_end", {})
        writer.write("run_end", {})


def test_transform(tmp_path):
    _write(tmp_path / "in.pst")
    seen = []

    def func(telescope_id, columns, metadata):
        seen.append((telescope_id, metadata))
        keep = columns["x"] % 2 == 0
        return {name: values[keep] * np.float32(2) for name, values in columns.items()}

    transform(
        tmp_path / "in.pst",
        tmp_path / "out.pst",
        func,
        stage_id=0x444C0001,
        frame_out="camera",
        efficiency_status={
            "qe": {
                "state": "applied",
                "source": "test",
                "wavelength_nm": [300.0, 400.0],
                "value": [0.1, 0.2],
            }
        },
    )
    records = list(obdeect_stream.Reader(str(tmp_path / "out.pst")))
    assert [r["type"] for r in records] == [
        "stream", "run", "efficiency_status", "event", "telescope", "photons", "photons", "blob",
        "telescope_end", "event_end", "run_end",
    ]  # fmt: skip
    assert records[0]["metadata"]["producer"] == "obdeect.stream.stage 0x444c0001"
    assert records[2]["metadata"]["factors"] == "qe"
    assert records[2]["metadata"]["qe.state"] == "applied"
    np.testing.assert_array_equal(records[2]["metadata"]["qe.value"], [0.1, 0.2])
    assert records[4]["metadata"]["frame"] == obdeect_stream.FRAMES["camera"]
    assert records[7]["data"] == b"\x01\x02"
    assert records[8]["metadata"] == {"photons": 52.0}
    for record, first, rows in ((records[5], 0, 3), (records[6], 3, 1)):
        assert record["frame"] == obdeect_stream.FRAMES["camera"]
        assert record["first_index"] == first
        assert len(record["columns"]["x"]) == rows
    np.testing.assert_array_equal(records[5]["columns"]["x"], [0, 4, 8])
    assert [t for t, _ in seen] == [4, 4]
    info = seen[0][1]
    assert info["run"] == {"run_number": 5}
    assert info["event"] == {"event_number": 1}
    assert info["telescope"]["telescope_id"] == 4
    assert info["frame"] == obdeect_stream.FRAMES["ground"]
    assert info["stage_id"] == 0x444C0001


def test_transform_keeps_existing_efficiency_status_first(tmp_path):
    _write(
        tmp_path / "in.pst",
        efficiency_status={"factors": "atmosphere", "atmosphere.state": "applied"},
    )
    transform(
        tmp_path / "in.pst",
        tmp_path / "out.pst",
        lambda t, c, m: c,
        stage_id=1,
        frame_out=2,
        efficiency_status={"qe": {"state": "pending", "source": "x"}},
    )
    efficiency_statuses = [
        r
        for r in obdeect_stream.Reader(str(tmp_path / "out.pst"))
        if r["type"] == "efficiency_status"
    ]
    assert [r["metadata"]["factors"] for r in efficiency_statuses] == ["atmosphere", "qe"]


def test_transform_rejects_reapplied_factor(tmp_path):
    _write(tmp_path / "in.pst", efficiency_status={"factors": "qe", "qe.state": "applied"})
    with pytest.raises(ValueError, match="already applied"):
        transform(
            tmp_path / "in.pst",
            tmp_path / "out.pst",
            lambda t, c, m: c,
            stage_id=1,
            frame_out="camera",
            efficiency_status={"qe": {"state": "applied", "source": "x"}},
        )


def test_transform_rejects_bad_arguments(tmp_path):
    _write(tmp_path / "in.pst")
    for kwargs in (
        {"stage_id": -1, "frame_out": "camera"},
        {"stage_id": 1, "frame_out": "moon"},
        {"stage_id": 1, "frame_out": 9},
    ):
        with pytest.raises(ValueError):
            transform(tmp_path / "in.pst", tmp_path / "o.pst", lambda t, c, m: c, **kwargs)
    with pytest.raises(TypeError, match="float32"):
        transform(
            tmp_path / "in.pst",
            tmp_path / "o.pst",
            lambda t, c, m: {k: v.astype(np.float64) for k, v in c.items()},
            stage_id=1,
            frame_out="camera",
        )


def test_efficiency_status_record_validation():
    with pytest.raises(ValueError, match="state"):
        efficiency_status_record({"qe": {"state": "maybe", "source": "x"}})
    with pytest.raises(ValueError, match="together"):
        efficiency_status_record({"qe": {"state": "applied", "source": "x", "value": [1.0]}})
    with pytest.raises(ValueError, match="length"):
        efficiency_status_record({
            "qe": {"state": "applied", "source": "x", "wavelength_nm": [1.0, 2.0], "value": [1.0]}
        })
    with pytest.raises(ValueError, match="unknown"):
        efficiency_status_record({"qe": {"state": "applied", "source": "x", "colour": 1}})


def test_conformance_streams():
    expected = json.loads((CONFORMANCE / "expected.json").read_text())
    records = list(obdeect_stream.Reader(str(CONFORMANCE / expected["valid"]["file"])))
    wanted = expected["valid"]["records"]
    assert len(records) == len(wanted)
    for record, want in zip(records, wanted, strict=True):
        assert record["type"] == want["type"]
        if "metadata" in want:
            assert set(record["metadata"]) == set(want["metadata"])
            for key, value in want["metadata"].items():
                np.testing.assert_array_equal(record["metadata"][key], value)
        elif want["type"] == "photons":
            assert record["telescope_id"] == want["telescope_id"]
            assert record["frame"] == want["frame"]
            assert record["first_index"] == want["first_index"]
            assert set(record["columns"]) == set(want["columns"])
            for name, values in want["columns"].items():
                np.testing.assert_array_equal(record["columns"][name], np.float32(values))
        else:
            assert record["eventio_type"] == want["eventio_type"]
            assert record["eventio_id"] == want["eventio_id"]
            assert record["data"] == bytes(want["data"])
    for name, message in expected["invalid"].items():
        with pytest.raises(Exception, match=re.escape(message)):
            list(obdeect_stream.Reader(str(CONFORMANCE / name)))


def test_transform_rejects_reapplied_factor_in_telescope_efficiency_status(tmp_path):
    with obdeect_stream.Writer(str(tmp_path / "in.pst"), "pytest") as writer:
        writer.write("run", {})
        writer.write("event", {})
        writer.write("telescope", {"telescope_id": 1, "frame": 1})
        writer.write("efficiency_status", {"factors": "qe", "qe.state": "applied"})
        writer.write("telescope_end", {})
        writer.write("event_end", {})
        writer.write("run_end", {})
    with pytest.raises(ValueError, match="already applied"):
        transform(
            tmp_path / "in.pst",
            tmp_path / "out.pst",
            lambda t, c, m: c,
            stage_id=1,
            frame_out="camera",
            efficiency_status={"qe": {"state": "applied", "source": "x"}},
        )


def test_transform_rejects_same_file_and_alias(tmp_path):
    source = tmp_path / "in.pst"
    _write(source)
    original = source.read_bytes()
    link = tmp_path / "link.pst"
    link.symlink_to(source)
    for target in (source, link):
        with pytest.raises(ValueError, match="same file"):
            transform(source, target, lambda tel, cols, info: cols, stage_id=1, frame_out="camera")
    assert source.read_bytes() == original


def test_failed_writer_has_no_stream_end(tmp_path):
    path = tmp_path / "partial.pst"
    with pytest.raises(RuntimeError, match="boom"):
        with obdeect_stream.Writer(str(path), "pytest") as writer:
            writer.write("run", {"run_number": 5})
            raise RuntimeError("boom")
    with pytest.raises(RuntimeError, match="truncated"):
        list(obdeect_stream.Reader(str(path)))


def test_transform_indices_continue_across_array_reuses(tmp_path):
    source = tmp_path / "reuse.pst"
    with obdeect_stream.Writer(str(source), "pytest") as writer:
        writer.write("run", {"run_number": 5})
        writer.write("event", {"event_number": 1})
        for array_index, first_index in ((0, 0), (1, 6)):
            writer.write(
                "telescope", {"telescope_id": 4, "frame": "ground", "array_index": array_index}
            )
            writer.write_photons(4, "ground", first_index, _chunk(6))
            writer.write("telescope_end", {"telescope_id": 4, "photons": 6.0})
        writer.write("event_end", {})
        writer.write("run_end", {})
    output = tmp_path / "out.pst"
    transform(source, output, lambda tel, cols, info: cols, stage_id=1, frame_out="ground")
    firsts = [
        r["first_index"] for r in obdeect_stream.Reader(str(output)) if r["type"] == "photons"
    ]
    assert firsts == [0, 6]
