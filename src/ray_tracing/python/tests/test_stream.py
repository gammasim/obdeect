"""Photon-stream stage: frames, focus, telescope filter, efficiency_status, determinism and records.

Tests that write or read photon streams are skipped when obdeect.stream is not installed.
"""

import math

import numpy as np
import pytest
from obdeect.ray_tracing import OpticalModel
from obdeect.ray_tracing.compiler import model_hash
from obdeect.ray_tracing.formats import ModelError
from obdeect.ray_tracing.model import _obscurer_support
from obdeect.ray_tracing.psf import containment
from obdeect.ray_tracing.stream import (
    CONSTANT_FACTOR,
    ground_to_telescope_rotation,
    lower_input_plane_height_m,
    minimum_input_plane_height_m,
    trace_stream,
)
from telescopes import dual, segmented, window_records

FOCAL_LENGTH_M = 8.0  # focal length of segmented()
RUN = {"run_number": 7, "wavelength_min_nm": 300.0, "wavelength_max_nm": 600.0}
SEED = 20261107


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
            assert height - np.dot(normal, model.axes_offset_m) > _obscurer_support(
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
    # Spot size: the radius holding 80 % of the photons is about 0.2 mm for this fixture (measured
    # with the segmented model at 4000 photons); 1 mm bounds it with a factor-of-5 margin.
    spot = containment(columns["x"], columns["y"], fraction=0.8)
    assert spot["radius_m"] < 1e-3
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


def test_mixed_directions_do_not_depend_on_chunking_or_threads(model, tmp_path, obdeect_stream):
    down = with_direction(beam(1500, (0.0, 0.0), radius_m=0.5), 0.0, 0.0, -1.0)
    up = with_direction(beam(1500, (0.0, 0.01), radius_m=0.5, seed=8), 0.0, 0.01, +1.0)
    photons = {name: np.concatenate([down[name], up[name]]) for name in down}
    outputs = []
    for index, (chunk, threads) in enumerate([(None, 1), (333, 2), (7, 1)]):
        source = tmp_path / f"m{index}.in"
        write_chunks(source, {1: photons}, chunk, obdeect_stream)
        output = tmp_path / f"m{index}.out"
        trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0), threads=threads)
        outputs.append(camera_photons(output, obdeect_stream))
    assert len(outputs[0]["x"]) > 0
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


def test_output_chunks_are_capped_and_numbered(model, tmp_path, obdeect_stream):
    # Bunches of 3 photons: more output rows than one chunk holds.
    photons = beam(200000, (0.0, 0.0), radius_m=6.0)
    photons["w"] = np.full(200000, 3.0, np.float32)
    source = write_chunks_of(tmp_path / "cap.in", photons, 50000, obdeect_stream)
    output = tmp_path / "cap.out"
    trace_stream(model, source, output, seed=SEED, pointing=(0, 0))
    chunks = [r for r in read_records(output, obdeect_stream) if r["type"] == "photons"]
    assert len(chunks) > 1
    expected_index = 0
    for chunk in chunks:
        assert len(chunk["columns"]["w"]) <= obdeect_stream.max_chunk_rows
        assert chunk["first_index"] == expected_index
        expected_index += len(chunk["columns"]["w"])
    telescope_end = [
        r for r in read_records(output, obdeect_stream) if r["type"] == "telescope_end"
    ]
    assert telescope_end[0]["metadata"]["photons"] == pytest.approx(
        float(camera_photons(output, obdeect_stream)["w"].sum()), rel=1e-6
    )


def test_output_index_continues_across_blocks_and_resets_per_event(model, tmp_path, obdeect_stream):
    """CORSIKA array reuse gives one event several blocks of the same telescope."""
    photons = beam(500, (0.0, 0.0))
    source = tmp_path / "reuse.in"
    with obdeect_stream.Writer(str(source), "pytest") as writer:
        writer.write("run", RUN)
        for event_number in (1, 2):
            writer.write("event", {"event_number": event_number})
            for first_index in (0, 500):
                writer.write("telescope", {"telescope_id": 1, "frame": 1})
                writer.write_photons(1, 1, first_index, photons)
                writer.write("telescope_end", {"telescope_id": 1, "photons": 500.0})
            writer.write("event_end", {})
        writer.write("run_end", {})
    output = tmp_path / "reuse.out"
    trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))
    events: list[list[dict]] = []
    for record in read_records(output, obdeect_stream):
        if record["type"] == "event":
            events.append([])
        elif record["type"] == "photons":
            events[-1].append(record)
    assert len(events) == 2
    for chunks in events:
        expected_index = 0
        for chunk in chunks:
            assert chunk["first_index"] == expected_index
            expected_index += len(chunk["columns"]["w"])
        assert expected_index > 0


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


@pytest.fixture
def window_model():
    records, files = window_records()
    return OpticalModel.from_records(
        records,
        files.__getitem__,
        configuration_seed=1,
        propagation_group_index=1.0,
        camera_window="sphere",
    )


def test_traced_camera_window_is_recorded_as_applied(model, window_model, tmp_path, obdeect_stream):
    _, _, plain_out = focus(
        model, tmp_path, obdeect_stream, beam(10, (0, 0)), pointing=(0, 0), name="nowin"
    )
    assert "camera_window" not in efficiency_status_of(plain_out, obdeect_stream)[0]["factors"]
    _, _, out = focus(
        window_model, tmp_path, obdeect_stream, beam(10, (0, 0)), pointing=(0, 0), name="win"
    )
    status = efficiency_status_of(out, obdeect_stream)[0]
    assert status["factors"].split(",") == [
        "mirror_reflectivity",
        "telescope_transmission",
        "camera_window",
    ]
    assert status["camera_window.state"] == "applied"
    assert status["camera_window.source"].endswith(window_model.sha256)


def test_camera_window_already_applied_is_rejected(window_model, tmp_path, obdeect_stream):
    applied = {"factors": "camera_window", "camera_window.state": "applied"}
    with pytest.raises(ModelError, match="camera_window is already applied"):
        focus(
            window_model,
            tmp_path,
            obdeect_stream,
            beam(10, (0, 0)),
            pointing=(0, 0),
            name="winapplied",
            run_efficiency_statuses=[applied],
        )


def test_angular_transmission_needs_a_model_with_transmission(model, tmp_path, obdeect_stream):
    records, files = segmented(telescope_transmission=None)
    bare = OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=1, propagation_group_index=1.0
    )
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
    with pytest.raises(ModelError, match="full or angular"):
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


def dual_without_obscurers():
    """SST-like dual mirror with no obscurers: its dish reaches above the camera at grazing."""
    records, files = dual()
    model = OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=1, propagation_group_index=1.0
    )
    document = model.document
    document["optics"]["obscurers"] = []
    document["sha256"] = model_hash(document)
    return OpticalModel(document)


def parallel_beam_in_both_frames(model, zenith_deg, azimuth_deg, count, radius_m, seed=5):
    """The same rays as telescope-frame and as ground-frame photons, for trace-stream.

    Rays move along the optical axis, in the telescope frame as (u, v, 0) moving along -z. The
    ground-frame photon of each ray is where it crosses the ground plane z = 0 (the telescope
    sphere centre plane), with the ground-frame direction. Both describe the same rays.
    """
    rng = np.random.default_rng(seed)
    radius = radius_m * np.sqrt(rng.random(count))
    angle = 2.0 * np.pi * rng.random(count)
    u, v = radius * np.cos(angle), radius * np.sin(angle)
    rotation = ground_to_telescope_rotation(zenith_deg, azimuth_deg)  # p_t = R p_g - O
    down = np.array([0.0, 0.0, -1.0])
    origin = np.stack([u, v, np.zeros(count)], axis=1) + np.asarray(model.axes_offset_m)
    # p_g = R^T (p_t + O); the ray crosses z_g = 0 where R[:, 2] . (origin + s down) = 0.
    s = -(origin @ rotation[:, 2]) / (down @ rotation[:, 2])
    ground_xy = ((origin + s[:, None] * down) @ rotation)[:, :2]
    ground_direction = down @ rotation
    zeros = np.zeros(count, np.float32)
    ones = np.ones(count, np.float32)
    wavelength = np.full(count, 400.0, np.float32)
    telescope = {
        "x": u.astype(np.float32),
        "y": v.astype(np.float32),
        "z": zeros,
        "cx": zeros,
        "cy": zeros,
        "t": zeros,
        "w": ones,
        "wavelength": wavelength,
    }
    ground = {
        "x": ground_xy[:, 0].astype(np.float32),
        "y": ground_xy[:, 1].astype(np.float32),
        "cx": np.full(count, ground_direction[0], np.float32),
        "cy": np.full(count, ground_direction[1], np.float32),
        "t": zeros,
        "w": ones,
        "wavelength": wavelength,
    }
    return telescope, ground


@pytest.mark.parametrize("zenith", [20.0, 60.0, 80.0, 86.775939])
def test_grazing_pointing_keeps_the_photons_of_the_same_rays(zenith, tmp_path, obdeect_stream):
    # At grazing pointing the dish reaches above the obscurer-based input plane. Photons must
    # start above the dish, not only above the obscurers, or the ones behind it never reach it.
    model = dual_without_obscurers()
    telescope, ground = parallel_beam_in_both_frames(model, zenith, 132.4, 20000, 2.0)
    reference, _, _ = focus(model, tmp_path, obdeect_stream, telescope, frame=2, name="t.pst")
    counts, _, _ = focus(
        model, tmp_path, obdeect_stream, ground, pointing=(zenith, 132.4), name="g.pst"
    )
    assert reference["photons_out"] > 0.5 * 20000
    assert counts["photons_in"] == 20000
    assert counts["photons_out"] == pytest.approx(reference["photons_out"], rel=1e-3)


def test_input_plane_clears_the_mirror_at_grazing_pointing():
    model = dual_without_obscurers()
    normal = ground_to_telescope_rotation(86.775939, 132.4)[:, 2]
    height = minimum_input_plane_height_m(model, 86.775939, 132.4)
    # The plane (ground height `height` above the sphere centre) lies above the mirror's
    # highest point along the ground normal, measured from the axes offset.
    support = dict(model.mirrors_support(normal))["primary"]
    assert height - float(np.dot(normal, model.axes_offset_m)) > support


@pytest.fixture(scope="module")
def dual_model():
    records, files = dual()
    return OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=1, propagation_group_index=1.0
    )


def two_telescope_stream(path):
    return write_stream(
        path,
        {1: beam(800, (0.0, 0.0), seed=1), 2: beam(600, (0.0, 0.0), radius_m=1.5, seed=2)},
        run_efficiency_statuses=[{"factors": "qe", "qe.state": "pending"}],
    )


def test_models_by_telescope_trace_as_separate_runs(model, dual_model, tmp_path, obdeect_stream):
    source = two_telescope_stream(tmp_path / "array.pst")
    output = tmp_path / "array.out.pst"

    counts = trace_stream(
        {1: model, 2: dual_model}, source, output, seed=SEED, pointing=(0.0, 0.0), threads=3
    )

    assert counts["photons_in"] == 1400
    for telescope, single in ((1, model), (2, dual_model)):
        alone = tmp_path / f"alone{telescope}.pst"
        trace_stream(single, source, alone, seed=SEED, pointing=(0.0, 0.0), telescopes=[telescope])
        expected = camera_photons(alone, obdeect_stream)
        traced = camera_photons(output, obdeect_stream, telescope)
        assert expected.keys() == traced.keys()
        for name, values in expected.items():
            assert np.array_equal(values, traced[name]), name


def test_models_by_telescope_mark_the_status_in_each_block(
    model, dual_model, tmp_path, obdeect_stream
):
    source = two_telescope_stream(tmp_path / "status.pst")
    output = tmp_path / "status.out.pst"

    trace_stream({1: model, 2: dual_model}, source, output, seed=SEED, pointing=(0.0, 0.0))

    records = read_records(output, obdeect_stream)
    first_event = [r["type"] for r in records].index("event")
    run_statuses = [
        r["metadata"] for r in records[:first_event] if r["type"] == "efficiency_status"
    ]
    assert run_statuses == [{"factors": "qe", "qe.state": "pending"}]
    block_statuses = {}
    for previous, record in zip(records, records[1:]):
        if previous["type"] == "telescope" and record["type"] == "efficiency_status":
            block_statuses[previous["metadata"]["telescope_id"]] = record["metadata"]
    assert block_statuses.keys() == {1, 2}
    for telescope, single in ((1, model), (2, dual_model)):
        status = block_statuses[telescope]
        assert status["mirror_reflectivity.state"] == "applied"
        assert status["mirror_reflectivity.source"].endswith(f"model {single.sha256}")


def test_models_by_telescope_drop_telescopes_without_a_model(model, tmp_path, obdeect_stream):
    source = two_telescope_stream(tmp_path / "drop.pst")
    output = tmp_path / "drop.out.pst"

    counts = trace_stream({2: model}, source, output, seed=SEED, pointing=(0.0, 0.0))

    records = read_records(output, obdeect_stream)
    assert [r["metadata"]["telescope_id"] for r in records if r["type"] == "telescope"] == [2]
    assert counts["photons_in"] == 600


def test_models_by_telescope_need_a_model_for_each_selected_telescope(model, tmp_path):
    source = two_telescope_stream(tmp_path / "missing.pst")
    with pytest.raises(ModelError, match=r"telescopes without an optical model: \[1\]"):
        trace_stream({2: model}, source, tmp_path / "missing.out.pst", seed=SEED, telescopes=[1, 2])
    with pytest.raises(ModelError, match="at least one optical model"):
        trace_stream({}, source, tmp_path / "missing.out.pst", seed=SEED)


def test_telescope_ids_outside_the_range_are_rejected(model, tmp_path):
    source = two_telescope_stream(tmp_path / "range.pst")
    output = tmp_path / "range.out.pst"
    for models in ({0: model}, {65536: model}):
        with pytest.raises(ModelError, match="outside 1..65535"):
            trace_stream(models, source, output, seed=SEED)
    with pytest.raises(ModelError, match="outside 1..65535"):
        trace_stream(model, source, output, seed=SEED, telescopes=[0])


def test_telescope_status_marks_own_factors_once_per_block(model, tmp_path, obdeect_stream):
    pending = {
        "factors": "qe,mirror_reflectivity",
        "qe.state": "pending",
        "mirror_reflectivity.state": "pending",
    }
    source = write_stream(
        tmp_path / "once.pst",
        {1: beam(200, (0.0, 0.0))},
        telescope_efficiency_statuses={1: pending},
    )
    output = tmp_path / "once.out.pst"
    trace_stream({1: model}, source, output, seed=SEED, pointing=(0, 0))

    records = read_records(output, obdeect_stream)
    types = [r["type"] for r in records]
    start = types.index("telescope")
    end = types.index("photons", start)
    statuses = [r["metadata"] for r in records[start:end] if r["type"] == "efficiency_status"]
    assert len(statuses) == 1
    names = statuses[0]["factors"].split(",")
    assert len(names) == len(set(names))
    assert statuses[0]["mirror_reflectivity.state"] == "applied"
    assert statuses[0]["qe.state"] == "pending"


def with_direction(photons, cx, cy, cz_sign):
    """Photons with the direction (cx, cy, cz_sign sqrt(1 - cx^2 - cy^2)) and a cz column."""
    count = len(photons["w"])
    cz = cz_sign * math.sqrt(1.0 - cx * cx - cy * cy)
    return {
        **photons,
        "cx": np.full(count, cx, np.float32),
        "cy": np.full(count, cy, np.float32),
        "cz": np.full(count, cz, np.float32),
    }


def test_photons_from_below_are_not_detected(model, tmp_path, obdeect_stream):
    # The dish is seen from its back: the photons are blocked, or reach the focal plane from below.
    photons = with_direction(beam(4000, (0.0, 0.0), radius_m=3.0), 0.0, 0.0, +1.0)
    source = write_stream(tmp_path / "up.in", {1: photons})
    counts = trace_stream(model, source, tmp_path / "up.out", seed=SEED, pointing=(0.0, 0.0))
    assert counts["photons_in"] == 4000
    assert counts["photons_out"] == 0


def test_mixed_directions_keep_the_downward_output(model, tmp_path, obdeect_stream):
    down = with_direction(beam(4000, (0.0, 0.0), radius_m=0.5), 0.0, 0.0, -1.0)
    up = with_direction(beam(2000, (0.01, 0.0), radius_m=0.5, seed=8), 0.01, 0.0, +1.0)
    # The upward photons follow the downward ones in the same chunk: their ids do not shift
    # the downward photons' random numbers, and the camera photons are the same.
    mixed = {name: np.concatenate([down[name], up[name]]) for name in down}
    alone = trace_stream(
        model,
        write_stream(tmp_path / "a.in", {1: down}),
        tmp_path / "a.out",
        seed=SEED,
        pointing=(0.0, 0.0),
    )
    both = trace_stream(
        model,
        write_stream(tmp_path / "b.in", {1: mixed}),
        tmp_path / "b.out",
        seed=SEED,
        pointing=(0.0, 0.0),
    )
    assert both["photons_in"] == 6000 and alone["photons_in"] == 4000
    assert both["photons_out"] == alone["photons_out"] > 0
    a = camera_photons(tmp_path / "a.out", obdeect_stream)
    b = camera_photons(tmp_path / "b.out", obdeect_stream)
    for name in a:
        assert a[name].tobytes() == b[name].tobytes(), name


def test_horizontal_photons_are_rejected(model, tmp_path, obdeect_stream):
    photons = with_direction(beam(10, (1.0, 0.0)), 1.0, 0.0, 1.0)
    photons["cz"] = np.zeros(10, np.float32)
    source = write_stream(tmp_path / "h.in", {1: photons})
    with pytest.raises(ModelError, match="horizontal photons"):
        trace_stream(model, source, tmp_path / "h.out", seed=SEED, pointing=(0.0, 0.0))


def test_lower_input_plane_is_below_every_optic(model):
    normal = np.array([0.0, 0.0, 1.0])
    bottoms = dict(model.optics_bottoms(normal))
    assert set(bottoms) >= {"primary"}
    # Pointing at the zenith: the plane lies a margin below the lowest optic (no axes offset here).
    plane = lower_input_plane_height_m(model, 0.0)
    assert plane == pytest.approx(min(bottoms.values()) - 1.0 + float(normal @ model.axes_offset_m))
    # A tilted pointing: the same plane along its zenith, below the optics' bottoms.
    zenith = ground_to_telescope_rotation(60.0, 0.0)[:, 2]
    height = lower_input_plane_height_m(model, 60.0)
    for name, bottom in model.optics_bottoms(zenith):
        assert height + float(np.dot(zenith, model.axes_offset_m)) < bottom, name


def test_upward_bunches_need_a_lower_plane_below_the_optics(model):
    bunches = {
        "x": np.zeros(1, np.float32),
        "y": np.zeros(1, np.float32),
        "cx": np.zeros(1, np.float32),
        "cy": np.zeros(1, np.float32),
        "cz": np.ones(1, np.float32),
        "time_ns": np.zeros(1, np.float32),
        "zem": np.zeros(1, np.float32),
        "size": np.ones(1, np.float32),
        "wavelength_nm": np.full(1, 400.0, np.float32),
        "bunch_id": np.arange(1, dtype=np.uint64),
    }
    common = {"frame_rotation": np.eye(3), "frame_offset_m": np.zeros(3), "length_unit_m": 1.0}
    tops = [z for z in (model.max_incoming_obscurer_z_m(), model.max_mirror_z_m()) if z is not None]
    plane = max(tops, default=0.0) + 1.0
    bottom = min((b for _, b in model.optics_bottoms(np.array([0.0, 0.0, 1.0]))), default=0.0)
    lower = bottom - 1.0
    model.trace_bunches(bunches, seed=1, plane_z=plane, lower_plane_z=lower, **common)
    with pytest.raises(ValueError, match="from below"):
        model.trace_bunches(bunches, seed=1, plane_z=plane, lower_plane_z=plane, **common)
    with pytest.raises(ValueError):
        model.trace_bunches(bunches, seed=1, plane_z=plane, **common)


def test_sources_from_below_emit_upward_photons(model):
    photons = model.star_photons(50, seed=3, towards_source=(0.0, 0.0, -1.0))
    assert np.all(photons["direction"][:, 2] > 0.0)
    # The photons start on a plane below every optic, the lowest one less the margin.
    bottom = min(b for _, b in model.optics_bottoms(np.array([0.0, 0.0, 1.0])))
    assert np.allclose(photons["position_m"][:, 2], bottom - 1.0, atol=1e-6)
    with pytest.raises(ValueError):
        model.star_photons(5, seed=3, towards_source=(1.0, 0.0, 0.0))


@pytest.mark.parametrize("event_number", [1, 255, 256, 300, 2**32 - 1])
def test_event_numbers_trace_with_their_own_random_numbers(
    model, tmp_path, obdeect_stream, event_number
):
    # Bunch ids are the photon indices: the event number is in the random key. Event numbers from
    # 256 on once overflowed the 40-bit photon id of the key and failed.
    photons = beam(3000, (0.0, 0.0), radius_m=6.0)
    source = write_stream(
        tmp_path / f"event{event_number}.in", {1: photons}, event={"event_number": event_number}
    )
    output = tmp_path / f"event{event_number}.out"
    counts = trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))
    assert counts["photons_out"] > 0
    assert camera_photons(output, obdeect_stream)["x"].size == counts["photons_out"]


def test_events_draw_different_random_numbers(model, tmp_path, obdeect_stream):
    photons = beam(3000, (0.0, 0.0), radius_m=6.0)
    x = {}
    for event_number in (255, 256, 300, 2**32 - 1):
        source = write_stream(
            tmp_path / f"e{event_number}.in", {1: photons}, event={"event_number": event_number}
        )
        output = tmp_path / f"e{event_number}.out"
        trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))
        x[event_number] = camera_photons(output, obdeect_stream)["x"]
    assert not np.array_equal(x[255], x[256])
    assert not np.array_equal(x[300], x[2**32 - 1])


def write_chunks_of(path, photons, chunk, obdeect_stream):
    """A stream of one telescope whose photon records hold `chunk` rows each."""
    rows = len(photons["w"])
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", RUN)
        writer.write("event", {"event_number": 1})
        writer.write("telescope", {"telescope_id": 1, "frame": 1})
        for start in range(0, rows, chunk):
            part = {k: v[start : start + chunk] for k, v in photons.items()}
            writer.write_photons(1, 1, start, part)
        writer.write("telescope_end", {"telescope_id": 1, "photons": float(rows)})
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def test_photon_sum_in_telescope_end_does_not_depend_on_chunking(model, tmp_path, obdeect_stream):
    # The photons metadata is the sum of the output weights, in output order.
    photons = beam(20000, (0.0, 0.0), radius_m=6.0)
    photons["w"] = np.linspace(0.5, 1.5, 20000, dtype=np.float32)
    totals = []
    for chunk in (1000, 7777, 20000):
        source = write_chunks_of(tmp_path / f"sum{chunk}.in", photons, chunk, obdeect_stream)
        output = tmp_path / f"sum{chunk}.out"
        trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))
        (end,) = [r for r in read_records(output, obdeect_stream) if r["type"] == "telescope_end"]
        totals.append(end["metadata"]["photons"])
    assert totals[0] > 0.0
    assert totals[0] == totals[1] == totals[2]


def test_cz_column_accepts_float32_near_horizontal_directions(model, tmp_path, obdeect_stream):
    # A unit vector in float32 with |cz| ~ 6e-4: an absolute check on cz rejected it. Upwards
    # photons are accepted too; a cz that does not match cx and cy is rejected.
    for cz, accepted in ((-0.000614837, True), (0.000614837, True), (-0.05, False)):
        photons = beam(1, (0.989117742, -0.147124648))
        photons["cz"] = np.array([cz], dtype=np.float32)
        source = write_stream(tmp_path / "cz.in", {1: photons})
        output = tmp_path / "cz.out"
        if not accepted:
            with pytest.raises(ModelError, match="cz"):
                trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))
            continue
        counts = trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))
        assert counts["photons_in"] == 1


def test_acceptance_column_is_rejected(model, tmp_path):
    photons = beam(10, (0.0, 0.0))
    photons["p"] = np.ones(10, np.float32)
    source = write_stream(tmp_path / "p.in", {1: photons})
    with pytest.raises(ModelError, match="photon column p"):
        trace_stream(model, source, tmp_path / "p.out", seed=SEED, pointing=(0.0, 0.0))


def test_single_wavelength_run_draws_from_one_nanometre_above_that_wavelength(
    model, tmp_path, obdeect_stream
):
    # One wavelength: the 1/lambda^2 law over [lo, lo + 1) nm (a wavelength-0 photon is drawn).
    photons = beam(300, (0.0, 0.0), radius_m=6.0)
    photons["wavelength"] = np.zeros(300, np.float32)
    run = dict(RUN, wavelength_min_nm=400.0, wavelength_max_nm=400.0)
    source = write_stream(tmp_path / "single.in", {1: photons}, run=run)
    output = tmp_path / "single.out"
    trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0))
    wavelength = camera_photons(output, obdeect_stream)["wavelength"]
    assert wavelength.size > 0
    assert np.all((wavelength >= 400.0) & (wavelength < 401.0))
