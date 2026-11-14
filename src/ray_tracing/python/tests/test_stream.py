"""Photon-stream stage: frames, focus, telescope filter, efficiency_status, determinism and records.

Tests that write or read photon streams are skipped when obdeect.stream is not installed.
"""

import math

import numpy as np
import pytest
from obdeect.ray_tracing import OpticalModel
from obdeect.ray_tracing import stream as stream_module
from obdeect.ray_tracing.cli import main
from obdeect.ray_tracing.formats import ModelError
from obdeect.ray_tracing.stream import (
    CONSTANT_FACTOR,
    ground_to_telescope_rotation,
    minimum_input_plane_height_m,
    trace_stream,
)
from telescopes import segmented

FOCAL_LENGTH_M = 8.0  # focal length of segmented()
RUN = {"run_number": 7, "wavelength_min_nm": 300.0, "wavelength_max_nm": 600.0}
SEED = 20261107


@pytest.fixture
def obdeect_stream():
    return pytest.importorskip("obdeect.stream")


@pytest.fixture(scope="module")
def model():
    records, files = segmented()
    return OpticalModel.from_records(records, files.__getitem__, configuration_seed=1)


def beam(count, direction, *, radius_m=0.5, wavelength_nm=400.0, seed=3):
    """Parallel photons on a disc in the ground plane z = 0, moving along `direction`."""
    rng = np.random.default_rng(seed)
    radius = radius_m * np.sqrt(rng.random(count))
    angle = 2.0 * np.pi * rng.random(count)
    return {
        "x": (radius * np.cos(angle)).astype(np.float32),
        "y": (radius * np.sin(angle)).astype(np.float32),
        "cx": np.full(count, direction[0], np.float32),
        "cy": np.full(count, direction[1], np.float32),
        "t": np.zeros(count, np.float32),
        "w": np.ones(count, np.float32),
        "wavelength": np.full(count, wavelength_nm, np.float32),
    }


def write_stream(
    path,
    photons,
    *,
    frame=1,
    run=RUN,
    event=None,
    run_efficiency_statuses=(),
    telescope_efficiency_statuses=None,
    blob=False,
):
    """Write a stream: `photons` maps telescope ids to photon columns (one chunk per call)."""
    if event is None:
        event = {"event_number": 1}
    obdeect_stream = pytest.importorskip("obdeect.stream")
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", run)
        for efficiency_status in run_efficiency_statuses:
            writer.write("efficiency_status", efficiency_status)
        if blob:
            writer.write_blob(1204, 7, b"\x01\x02\x03")
        writer.write("event", event)
        for telescope_id, columns in photons.items():
            writer.write("telescope", {"telescope_id": telescope_id, "frame": frame})
            if telescope_efficiency_statuses and telescope_id in telescope_efficiency_statuses:
                writer.write("efficiency_status", telescope_efficiency_statuses[telescope_id])
            writer.write_photons(telescope_id, frame, 0, columns)
            writer.write(
                "telescope_end",
                {"telescope_id": telescope_id, "photons": float(columns["w"].sum())},
            )
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def read_records(path, obdeect_stream):
    return list(obdeect_stream.Reader(str(path)))


def camera_photons(path, obdeect_stream, telescope_id=None):
    """Concatenated C-frame columns of the photon records of one telescope (or all)."""
    chunks = [
        r
        for r in read_records(path, obdeect_stream)
        if r["type"] == "photons" and telescope_id in (None, r["telescope_id"])
    ]
    names = chunks[0]["columns"].keys()
    return {name: np.concatenate([c["columns"][name] for c in chunks]) for name in names}


def efficiency_status_of(path, obdeect_stream):
    return [
        r["metadata"]
        for r in read_records(path, obdeect_stream)
        if r["type"] == "efficiency_status"
    ]


def focus(
    model,
    tmp_path,
    obdeect_stream,
    photons,
    *,
    pointing=None,
    name="s.pst",
    telescope_transmission="full",
    **options,
):
    source = write_stream(tmp_path / f"{name}.in", {1: photons}, **options)
    output = tmp_path / f"{name}.out"
    counts = trace_stream(
        model,
        source,
        output,
        seed=SEED,
        pointing=pointing,
        telescope_transmission=telescope_transmission,
    )
    columns = camera_photons(output, obdeect_stream)
    return counts, columns, output


def test_rotation_maps_the_optical_axis_to_minus_z():
    theta, phi = math.radians(12.0), math.radians(73.0)
    # CORSIKA propagation direction of a shower along the optical axis.
    direction = np.array([
        -math.sin(theta) * math.cos(phi),
        math.sin(theta) * math.sin(phi),
        -math.cos(theta),
    ])
    rotation = ground_to_telescope_rotation(12.0, 73.0)
    assert np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-12)
    assert np.allclose(rotation @ direction, [0.0, 0.0, -1.0], atol=1e-12)


def test_input_plane_is_above_the_incoming_obscurers(model):
    top = model.max_incoming_obscurer_z_m()
    assert minimum_input_plane_height_m(model, 0.0) == pytest.approx(top + 1.0)
    normal = ground_to_telescope_rotation(60.0, 0.0)[:, 2]
    height = minimum_input_plane_height_m(model, 60.0)
    for obscurer in model.document["optics"]["obscurers"]:
        if "incoming" in obscurer.get("legs", []):
            assert height - np.dot(normal, model.axes_offset_m) > stream_module._obscurer_support(
                obscurer, normal
            )
    with pytest.raises(ModelError):
        minimum_input_plane_height_m(model, 90.0)


def test_output_telescope_frame_is_camera(model, tmp_path, obdeect_stream):
    _, _, output = focus(model, tmp_path, obdeect_stream, beam(10, (0, 0)), pointing=(0, 0))
    telescope = next(r for r in read_records(output, obdeect_stream) if r["type"] == "telescope")
    assert telescope["metadata"]["frame"] == 3


def test_on_axis_beam_focuses_at_the_centre(model, tmp_path, obdeect_stream):
    counts, columns, _ = focus(
        model, tmp_path, obdeect_stream, beam(4000, (0.0, 0.0, -1.0)[:2]), pointing=(0.0, 0.0)
    )
    assert counts["photons_in"] == 4000
    assert 0 < counts["photons_out"] < 4000
    assert abs(columns["x"].mean()) < 1e-3
    assert abs(columns["y"].mean()) < 1e-3
    # Reflected photons travel from the dish up towards the camera: +z in the telescope frame.
    assert np.all(columns["cz"] > 0.0)
    assert np.allclose(columns["cx"] ** 2 + columns["cy"] ** 2 + columns["cz"] ** 2, 1.0, atol=1e-5)


def test_off_axis_beam_shifts_by_f_tan_theta(model, tmp_path, obdeect_stream):
    angle = math.radians(2.0)
    _, columns, _ = focus(
        model, tmp_path, obdeect_stream, beam(4000, (-math.sin(angle), 0.0)), pointing=(0.0, 0.0)
    )
    expected = FOCAL_LENGTH_M * math.tan(angle)
    assert abs(abs(columns["x"].mean()) - expected) < 0.02 * expected
    assert abs(columns["y"].mean()) < 1e-3


def test_pointing_rotates_the_ground_frame(model, tmp_path, obdeect_stream):
    # A beam parallel to a telescope tilted by 2 deg is focused at the centre.
    angle = math.radians(2.0)
    _, columns, _ = focus(
        model,
        tmp_path,
        obdeect_stream,
        beam(4000, (-math.sin(angle), 0.0)),
        pointing=(2.0, 0.0),
        name="tilted",
    )
    assert abs(columns["x"].mean()) < 1e-3
    # A vertical beam is then off axis by f tan(2 deg).
    _, columns, _ = focus(
        model,
        tmp_path,
        obdeect_stream,
        beam(4000, (0.0, 0.0)),
        pointing=(2.0, 0.0),
        name="vertical",
    )
    assert abs(abs(columns["x"].mean()) - FOCAL_LENGTH_M * math.tan(angle)) < 0.02 * (
        FOCAL_LENGTH_M * math.tan(angle)
    )


def test_pointing_is_read_from_the_run_record(model, tmp_path, obdeect_stream):
    photons = beam(2000, (0.0, 0.0))
    run = dict(RUN, telescope_pointing_zenith_rad=math.radians(2.0))
    run["telescope_pointing_azimuth_rad"] = math.radians(30.0)
    _, from_run, _ = focus(model, tmp_path, obdeect_stream, photons, run=run, name="run")
    _, from_args, _ = focus(
        model, tmp_path, obdeect_stream, photons, pointing=(2.0, 30.0), name="args"
    )
    assert np.array_equal(from_run["x"], from_args["x"])


def test_run_pointing_needs_both_keys(model, tmp_path, obdeect_stream):
    run = dict(RUN, telescope_pointing_zenith_rad=0.1)
    source = write_stream(tmp_path / "half.pst", {1: beam(10, (0.0, 0.0))}, run=run)
    with pytest.raises(ModelError, match="telescope_pointing_azimuth_rad"):
        trace_stream(model, source, tmp_path / "out.pst", seed=SEED)


@pytest.mark.parametrize("alias", ["same", "symlink", "hardlink"])
def test_trace_stream_rejects_input_aliases(model, tmp_path, obdeect_stream, alias):
    source = write_stream(tmp_path / "input.pst", {1: beam(10, (0.0, 0.0))})
    if alias == "same":
        output = source
    elif alias == "symlink":
        output = tmp_path / "alias.pst"
        output.symlink_to(source)
    else:
        output = tmp_path / "alias.pst"
        output.hardlink_to(source)

    with pytest.raises(ModelError, match="same file"):
        trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))


def test_ground_photons_need_a_pointing(model, tmp_path, obdeect_stream):
    # The event header carries no telescope pointing: its CORSIKA azimuth is ignored.
    event = {"event_number": 1, "zenith_rad": 0.3, "azimuth_rad": 1.0}
    source = write_stream(tmp_path / "nopointing.pst", {1: beam(10, (0.0, 0.0))}, event=event)
    with pytest.raises(ModelError, match="no telescope pointing"):
        trace_stream(model, source, tmp_path / "out.pst", seed=SEED)


def test_on_axis_beam_at_a_tilted_pointing_focuses_at_the_centre(model, tmp_path, obdeect_stream):
    # Shower from the optical axis of a telescope at zenith 20 deg, azimuth 0 (array frame):
    # propagation direction (-sin(zenith) cos(azimuth), sin(zenith) sin(azimuth), -cos(zenith)).
    direction = (-math.sin(math.radians(20.0)), 0.0)
    _, columns, _ = focus(
        model, tmp_path, obdeect_stream, beam(4000, direction), pointing=(20.0, 0.0), name="on"
    )
    assert len(columns["x"]) > 1000
    assert abs(columns["x"].mean()) < 1e-3
    assert abs(columns["y"].mean()) < 1e-3


def test_wrong_azimuth_convention_misses_the_camera(model, tmp_path, obdeect_stream):
    # Regression: the same beam with the CORSIKA propagation azimuth (about 175.68 deg for
    # an OBS-North shower at zenith 20 deg) as telescope azimuth is about 40 deg off axis.
    direction = (-math.sin(math.radians(20.0)), 0.0)
    _, columns, _ = focus(
        model, tmp_path, obdeect_stream, beam(4000, direction), pointing=(20.0, 175.68), name="bad"
    )
    assert math.hypot(columns["x"].mean(), columns["y"].mean()) > 1.0


def test_telescope_frame_matches_ground_frame_at_zenith(model, tmp_path, obdeect_stream):
    photons = beam(2000, (0.0, 0.0))
    _, ground, _ = focus(model, tmp_path, obdeect_stream, photons, pointing=(0.0, 0.0), name="g")
    telescope = dict(photons, z=np.zeros(2000, np.float32))
    _, frame_t, _ = focus(
        model, tmp_path, obdeect_stream, telescope, frame=2, pointing=None, name="t"
    )
    for name in ("x", "y", "cx", "cy", "cz", "t", "w", "wavelength"):
        assert np.array_equal(ground[name], frame_t[name]), name


def test_telescope_frame_needs_z(model, tmp_path, obdeect_stream):
    with pytest.raises(ModelError, match="column z"):
        focus(model, tmp_path, obdeect_stream, beam(10, (0.0, 0.0)), frame=2, name="noz")


def test_unsupported_frame_and_column_are_rejected(model, tmp_path, obdeect_stream):
    photons = dict(beam(10, (0.0, 0.0)), p=np.ones(10, np.float32))
    with pytest.raises(ModelError, match="column p"):
        focus(model, tmp_path, obdeect_stream, photons, pointing=(0.0, 0.0), name="p")
    with pytest.raises(ModelError, match="frame 0"):
        write_stream(tmp_path / "sky.pst", {1: beam(10, (0.0, 0.0))}, frame=0)
        trace_stream(model, tmp_path / "sky.pst", tmp_path / "x.pst", seed=SEED, pointing=(0, 0))


def test_unknown_photon_column_is_rejected():
    tracer = object.__new__(stream_module._Tracer)
    tracer.telescope = stream_module._Telescope(1, True)
    tracer.event = {"event_number": 1}
    with pytest.raises(ModelError, match="unsupported photon columns: extra"):
        tracer._photons({"telescope_id": 1, "columns": {"w": np.ones(1), "extra": np.ones(1)}})


def test_telescope_filter_drops_other_telescopes(model, tmp_path, obdeect_stream):
    source = write_stream(
        tmp_path / "two.pst",
        {1: beam(500, (0.0, 0.0), seed=1), 2: beam(300, (0.0, 0.0), seed=2)},
        blob=True,
        telescope_efficiency_statuses={1: {"factors": "qe", "qe.state": "pending"}},
    )
    output = tmp_path / "two.out.pst"
    counts = trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0), telescopes=[2])
    records = read_records(output, obdeect_stream)
    telescopes = [r["metadata"]["telescope_id"] for r in records if r["type"] == "telescope"]
    assert telescopes == [2]
    assert {r["telescope_id"] for r in records if r["type"] == "photons"} == {2}
    assert [r["metadata"]["telescope_id"] for r in records if r["type"] == "telescope_end"] == [2]
    assert not any(
        r["type"] == "efficiency_status" and r["metadata"].get("qe.state") for r in records[3:]
    )
    assert counts["photons_in"] == 300


def test_efficiency_status_pending_and_absent_apply_the_factors_in_full(
    model, tmp_path, obdeect_stream
):
    photons = beam(3000, (0.0, 0.0))
    _, absent, _ = focus(model, tmp_path, obdeect_stream, photons, pointing=(0, 0), name="absent")
    pending = {
        "factors": "mirror_reflectivity,telescope_transmission,qe",
        "mirror_reflectivity.state": "pending",
        "telescope_transmission.state": "pending",
        "qe.state": "pending",
    }
    _, with_pending, _ = focus(
        model,
        tmp_path,
        obdeect_stream,
        photons,
        pointing=(0, 0),
        name="pend",
        run_efficiency_statuses=[pending],
    )
    assert np.array_equal(absent["w"], with_pending["w"])
    assert np.array_equal(absent["x"], with_pending["x"])


def test_efficiency_status_envelope_divides_by_the_envelope(model, tmp_path, obdeect_stream):
    photons = beam(3000, (0.0, 0.0))
    _, plain, _ = focus(model, tmp_path, obdeect_stream, photons, pointing=(0, 0), name="plain")
    # Envelope 2 on both factors at every wavelength; the traced photons carry 1/2.
    flat = {
        "factors": "mirror_reflectivity",
        "mirror_reflectivity.state": "envelope",
        "mirror_reflectivity.wavelength_nm": np.array([300.0, 600.0]),
        "mirror_reflectivity.value": np.array([2.0, 2.0]),
    }
    _, flat_out, _ = focus(
        model,
        tmp_path,
        obdeect_stream,
        photons,
        pointing=(0, 0),
        name="flat",
        run_efficiency_statuses=[flat],
    )
    # The same envelope written as a steeper table: at 400 nm both give 1/2.
    sloped = dict(flat, **{"mirror_reflectivity.value": np.array([1.0, 4.0])})
    _, sloped_out, _ = focus(
        model,
        tmp_path,
        obdeect_stream,
        photons,
        pointing=(0, 0),
        name="sloped",
        run_efficiency_statuses=[sloped],
    )
    assert np.array_equal(flat_out["x"], sloped_out["x"])
    assert np.array_equal(flat_out["w"], sloped_out["w"])
    assert np.all(flat_out["w"] <= 0.5 + 1e-6)
    # Unbiased: the weight sum matches the full application (roulette) within statistics.
    assert float(flat_out["w"].sum()) == pytest.approx(0.5 * plain["w"].sum(), rel=0.05)


def test_efficiency_status_applied_is_rejected(model, tmp_path, obdeect_stream):
    applied = {"factors": "telescope_transmission", "telescope_transmission.state": "applied"}
    with pytest.raises(ModelError, match="already applied"):
        focus(
            model,
            tmp_path,
            obdeect_stream,
            beam(10, (0, 0)),
            pointing=(0, 0),
            name="applied",
            run_efficiency_statuses=[applied],
        )


def test_telescope_efficiency_status_overrides_the_run_efficiency_status(
    model, tmp_path, obdeect_stream
):
    run_pending = {"factors": "mirror_reflectivity", "mirror_reflectivity.state": "pending"}
    tel_envelope = {
        "factors": "mirror_reflectivity",
        "mirror_reflectivity.state": "envelope",
        "mirror_reflectivity.wavelength_nm": np.array([300.0, 600.0]),
        "mirror_reflectivity.value": np.array([2.0, 2.0]),
    }
    source = write_stream(
        tmp_path / "tel.pst",
        {1: beam(2000, (0.0, 0.0))},
        run_efficiency_statuses=[run_pending],
        telescope_efficiency_statuses={1: tel_envelope},
    )
    output = tmp_path / "tel.out.pst"
    trace_stream(model, source, output, seed=SEED, pointing=(0, 0))
    assert np.all(camera_photons(output, obdeect_stream)["w"] <= 0.5 + 1e-6)


def test_efficiency_status_output_marks_own_factors_and_keeps_others(
    model, tmp_path, obdeect_stream
):
    qe = {
        "factors": "qe,mirror_reflectivity",
        "qe.state": "envelope",
        "qe.wavelength_nm": np.array([300.0, 600.0]),
        "qe.value": np.array([0.3, 0.4]),
        "qe.source": "sim_telarray",
        "mirror_reflectivity.state": "pending",
    }
    _, _, output = focus(
        model,
        tmp_path,
        obdeect_stream,
        beam(10, (0, 0)),
        pointing=(0, 0),
        name="out",
        run_efficiency_statuses=[qe],
    )
    records = read_records(output, obdeect_stream)
    types = [r["type"] for r in records]
    run_index = types.index("run")
    assert types[run_index + 1] == "efficiency_status"
    rewritten = records[run_index + 1]["metadata"]
    own = records[run_index + 2]["metadata"]
    assert own["mirror_reflectivity.state"] == "applied"
    assert own["telescope_transmission.state"] == "applied"
    assert rewritten["mirror_reflectivity.state"] == "pending"
    assert rewritten["qe.state"] == "envelope"
    assert np.array_equal(rewritten["qe.value"], np.array([0.3, 0.4]))
    assert rewritten["qe.source"] == "sim_telarray"


def test_invalid_envelope_table_is_rejected(model, tmp_path, obdeect_stream):
    bad = {
        "factors": "mirror_reflectivity",
        "mirror_reflectivity.state": "envelope",
        "mirror_reflectivity.wavelength_nm": np.array([600.0, 300.0]),
        "mirror_reflectivity.value": np.array([2.0, 2.0]),
    }
    with pytest.raises(ModelError, match="increasing"):
        focus(
            model,
            tmp_path,
            obdeect_stream,
            beam(10, (0, 0)),
            pointing=(0, 0),
            name="bad",
            run_efficiency_statuses=[bad],
        )


def test_output_does_not_depend_on_chunking_or_threads(model, tmp_path, obdeect_stream):
    photons = beam(2000, (-0.01, 0.005))
    outputs = []
    for index, (chunk, threads) in enumerate([(None, 1), (333, 2), (7, 1)]):
        source = tmp_path / f"c{index}.in"
        write_chunks(source, {1: photons}, chunk, obdeect_stream)
        output = tmp_path / f"c{index}.out"
        trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0), threads=threads)
        outputs.append(camera_photons(output, obdeect_stream))
    for other in outputs[1:]:
        for name in outputs[0]:
            assert np.array_equal(outputs[0][name], other[name]), name


def write_chunks(path, photons, chunk, obdeect_stream):
    """A stream whose photon chunks hold `chunk` rows each (None: one chunk)."""
    columns = next(iter(photons.values()))
    rows = len(columns["w"])
    step = chunk or rows
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", RUN)
        writer.write("event", {"event_number": 1, "zenith_rad": 0.0, "azimuth_rad": 0.0})
        for telescope_id, cols in photons.items():
            writer.write("telescope", {"telescope_id": telescope_id, "frame": 1})
            for start in range(0, rows, step):
                part = {k: v[start : start + step] for k, v in cols.items()}
                writer.write_photons(telescope_id, 1, start, part)
            writer.write("telescope_end", {"telescope_id": telescope_id, "photons": float(rows)})
        writer.write("event_end", {})
        writer.write("run_end", {})


def test_output_chunks_are_capped_and_numbered(model, tmp_path, obdeect_stream, monkeypatch):
    monkeypatch.setattr(stream_module, "MAX_CHUNK_ROWS", 100)
    _, _, output = focus(
        model, tmp_path, obdeect_stream, beam(2000, (0.0, 0.0)), pointing=(0, 0), name="cap"
    )
    chunks = [r for r in read_records(output, obdeect_stream) if r["type"] == "photons"]
    assert len(chunks) > 1
    expected_index = 0
    for chunk in chunks:
        assert len(chunk["columns"]["w"]) <= 100
        assert chunk["first_index"] == expected_index
        expected_index += len(chunk["columns"]["w"])
    telescope_end = [
        r for r in read_records(output, obdeect_stream) if r["type"] == "telescope_end"
    ]
    assert telescope_end[0]["metadata"]["photons"] == pytest.approx(
        float(camera_photons(output, obdeect_stream)["w"].sum()), rel=1e-6
    )


def test_records_pass_through(model, tmp_path, obdeect_stream):
    source = write_stream(
        tmp_path / "pass.pst",
        {1: beam(20, (0.0, 0.0))},
        blob=True,
        event={"event_number": 9, "zenith_rad": 0.0, "azimuth_rad": 0.0, "energy_gev": 1.5},
    )
    output = tmp_path / "pass.out.pst"
    trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))
    records = read_records(output, obdeect_stream)
    types = [r["type"] for r in records]
    assert types[0] == "stream"
    assert types[-1] == "run_end"
    assert types.count("blob") == 1
    assert "pytest" in records[0]["metadata"]["producer"]
    assert any(r["type"] == "blob" and r["data"] == b"\x01\x02\x03" for r in records)
    run = next(r for r in records if r["type"] == "run")["metadata"]
    assert run == RUN
    event = next(r for r in records if r["type"] == "event")["metadata"]
    assert event["energy_gev"] == 1.5


def test_cli_trace_stream(model, tmp_path, obdeect_stream, capsys):
    model_path = tmp_path / "model.json"
    model.save(model_path)
    source = write_stream(tmp_path / "cli.in.pst", {1: beam(500, (0.0, 0.0))})
    output = tmp_path / "cli.out.pst"
    code = main([
        "trace-stream",
        "--model",
        str(model_path),
        "--seed",
        "3",
        "--pointing-zenith-deg",
        "0",
        "--pointing-azimuth-deg",
        "0",
        "-i",
        str(source),
        "-o",
        str(output),
    ])
    assert code == 0
    assert "trace-stream:" in capsys.readouterr().err
    assert any(r["type"] == "photons" for r in read_records(output, obdeect_stream))


def test_cli_needs_a_pointing_for_ground_photons(model, tmp_path, capsys):
    model_path = tmp_path / "model.json"
    model.save(model_path)
    source = write_stream(tmp_path / "cli.nop.pst", {1: beam(10, (0.0, 0.0))})
    code = main([
        "trace-stream",
        "--model",
        str(model_path),
        "--seed",
        "3",
        "-i",
        str(source),
        "-o",
        str(tmp_path / "cli.nop.out.pst"),
    ])
    assert code == 1
    assert "--pointing-zenith-deg" in capsys.readouterr().err


def test_cli_rejects_half_a_pointing(model, tmp_path, capsys):
    model_path = tmp_path / "model.json"
    model.save(model_path)
    code = main([
        "trace-stream",
        "--model",
        str(model_path),
        "--seed",
        "3",
        "--pointing-zenith-deg",
        "5",
        "-i",
        str(tmp_path / "none.pst"),
    ])
    assert code == 1
    assert "pointing" in capsys.readouterr().err


def test_angular_transmission_drops_the_on_axis_value_from_on_axis_photons(
    model, tmp_path, obdeect_stream
):
    # Telescope transmission on_axis 0.9 (segmented). Photons on the axis see t = t0 in the full
    # trace and t / t0 = 1 in the angular one, so the expected sum of output weights (the count
    # of survivors, roulette on) grows by 1 / t0. The check is statistical (Bernoulli variance).
    on_axis = 0.9
    photons = beam(20000, (0.0, 0.0))
    _, full, _ = focus(model, tmp_path, obdeect_stream, photons, pointing=(0, 0), name="tfull")
    _, angular, _ = focus(
        model,
        tmp_path,
        obdeect_stream,
        photons,
        pointing=(0, 0),
        name="tang",
        telescope_transmission="angular",
    )
    sum_full = float(full["w"].sum(dtype=np.float64))
    sum_angular = float(angular["w"].sum(dtype=np.float64))
    variance = float(np.sum(full["w"].astype(np.float64) ** 2))
    variance += on_axis**2 * float(np.sum(angular["w"].astype(np.float64) ** 2))
    print(f"on-axis weight sums: {on_axis * sum_angular:.1f} vs {sum_full:.1f}")
    assert sum_full > 1000
    assert abs(on_axis * sum_angular - sum_full) < 5 * math.sqrt(variance)
    # The change is visible: the angular sum is about 1 / t0 larger.
    assert sum_angular > 1.05 * sum_full


def test_efficiency_status_records_the_transmission_part(model, tmp_path, obdeect_stream):
    _, _, full_out = focus(
        model, tmp_path, obdeect_stream, beam(10, (0, 0)), pointing=(0, 0), name="lfull"
    )
    full = efficiency_status_of(full_out, obdeect_stream)[0]
    assert full["telescope_transmission.part"] == "full"
    assert full["telescope_transmission.state"] == "applied"
    assert CONSTANT_FACTOR not in full["factors"]
    _, _, angular_out = focus(
        model,
        tmp_path,
        obdeect_stream,
        beam(10, (0, 0)),
        pointing=(0, 0),
        name="lang",
        telescope_transmission="angular",
    )
    angular = efficiency_status_of(angular_out, obdeect_stream)[0]
    assert angular["telescope_transmission.part"] == "angular"
    assert angular["telescope_transmission.state"] == "applied"
    assert angular["factors"] == (
        "mirror_reflectivity,telescope_transmission,telescope_transmission_constant"
    )
    assert angular[f"{CONSTANT_FACTOR}.state"] == "pending"
    assert float(angular[f"{CONSTANT_FACTOR}.value"]) == pytest.approx(0.9)


def test_angular_efficiency_status_marks_telescope_level_transmission(
    model, tmp_path, obdeect_stream
):
    source = write_stream(
        tmp_path / "tl.pst",
        {1: beam(10, (0.0, 0.0))},
        telescope_efficiency_statuses={
            1: {"factors": "telescope_transmission", "telescope_transmission.state": "pending"}
        },
    )
    output = tmp_path / "tl.out.pst"
    trace_stream(
        model,
        source,
        output,
        seed=SEED,
        pointing=(0, 0),
        telescope_transmission="angular",
    )
    run_efficiency_status, telescope_efficiency_status = efficiency_status_of(
        output, obdeect_stream
    )[:2]
    assert run_efficiency_status["telescope_transmission.part"] == "angular"
    assert telescope_efficiency_status["telescope_transmission.part"] == "angular"
    assert telescope_efficiency_status[f"{CONSTANT_FACTOR}.state"] == "pending"


def test_angular_transmission_needs_a_model_with_transmission(model, tmp_path, obdeect_stream):
    records, files = segmented(telescope_transmission=None)
    bare = OpticalModel.from_records(records, files.__getitem__, configuration_seed=1)
    with pytest.raises(ModelError, match="has none"):
        focus(
            bare,
            tmp_path,
            obdeect_stream,
            beam(10, (0, 0)),
            pointing=(0, 0),
            name="bare",
            telescope_transmission="angular",
        )
    # The full transmission has nothing to apply here and still runs.
    _, _, output = focus(
        bare, tmp_path, obdeect_stream, beam(10, (0, 0)), pointing=(0, 0), name="bare-full"
    )
    if not bare.has_structure:
        assert (
            "telescope_transmission"
            not in efficiency_status_of(output, obdeect_stream)[0]["factors"]
        )


def test_unknown_transmission_part_and_existing_constant_are_rejected(
    model, tmp_path, obdeect_stream
):
    with pytest.raises(ModelError, match="must be one of"):
        focus(
            model,
            tmp_path,
            obdeect_stream,
            beam(10, (0, 0)),
            pointing=(0, 0),
            name="unknown",
            telescope_transmission="partial",
        )
    named = {
        "factors": "telescope_transmission,telescope_transmission_constant",
        "telescope_transmission.state": "pending",
        f"{CONSTANT_FACTOR}.state": "pending",
        f"{CONSTANT_FACTOR}.value": 0.5,
    }
    constant_only = {"factors": CONSTANT_FACTOR, f"{CONSTANT_FACTOR}.state": "pending"}
    with pytest.raises(ModelError, match="already names"):
        focus(
            model,
            tmp_path,
            obdeect_stream,
            beam(10, (0, 0)),
            pointing=(0, 0),
            run_efficiency_statuses=[constant_only],
            telescope_transmission="angular",
            name="constant-only",
        )
    with pytest.raises(ModelError, match="already names"):
        focus(
            model,
            tmp_path,
            obdeect_stream,
            beam(10, (0, 0)),
            pointing=(0, 0),
            name="named",
            run_efficiency_statuses=[named],
            telescope_transmission="angular",
        )


def test_cli_telescope_transmission_option(model, tmp_path, obdeect_stream, capsys):
    model_path = tmp_path / "model.json"
    model.save(model_path)
    source = write_stream(tmp_path / "opt.in.pst", {1: beam(50, (0.0, 0.0))})
    base = [
        "trace-stream",
        "--model",
        str(model_path),
        "--seed",
        "3",
        "--pointing-zenith-deg",
        "0",
        "--pointing-azimuth-deg",
        "0",
        "-i",
        str(source),
        "-o",
    ]
    parts = {}
    for label, extra in (("default", []), ("angular", ["--telescope-transmission", "angular"])):
        output = tmp_path / f"opt.{label}.pst"
        assert main([*base, str(output), *extra]) == 0
        parts[label] = efficiency_status_of(output, obdeect_stream)[0][
            "telescope_transmission.part"
        ]
    assert parts == {"default": "full", "angular": "angular"}
    capsys.readouterr()
    with pytest.raises(SystemExit) as error:
        main([*base, str(tmp_path / "bad.pst"), "--telescope-transmission", "partial"])
    assert error.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


def test_single_wavelength_run_keys_give_that_wavelength():
    law = stream_module._wavelength_law({"wavelength_min_nm": 400.0, "wavelength_max_nm": 400.0})
    assert law == (400.0, 401.0)
    # Any draw of the 1/lambda^2 law over [400, 401) truncates to 400 nm.
    lo, hi = law
    for u in (1e-9, 0.5, 1 - 1e-9):
        assert math.floor(1.0 / (1.0 / lo - u * (1.0 / lo - 1.0 / hi))) == 400
    assert stream_module._wavelength_law({"wavelength_min_nm": 300, "wavelength_max_nm": 550}) == (
        300.0,
        550.0,
    )
    assert stream_module._wavelength_law({}) is None
