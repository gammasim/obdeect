"""Sphere runs (ground_representation sphere): G rows start at their own points.

A sphere row is the point where the photon's path enters the telescope sphere (or its emission
point inside it). The tracer starts the photon there: sphere rows of a plane run's photons give the
same detected photons; horizontal rows are allowed; the sphere must enclose the optical model.
Tests that write or read photon streams are skipped when obdeect.stream is not installed.
"""

import subprocess

import numpy as np
import pytest
from obdeect.ray_tracing import OpticalModel
from obdeect.ray_tracing.compiler import model_hash
from obdeect.ray_tracing.formats import ModelError
from obdeect.ray_tracing.stream import (
    enclosure_distances,
    ground_to_telescope_rotation,
    minimum_input_plane_height_m,
    trace_stream,
)
from telescopes import segmented
from test_stream import (
    RUN,
    SEED,
    camera_photons,
    write_stream,
)
from test_stream_native import _assert_same_stream, _native_argv, _native_program

SPEED_OF_LIGHT_M_PER_NS = 0.299792458
SPHERE_RADIUS_M = 40.0
SPHERE_RUN = {
    **RUN,
    "ground_representation": "sphere",
    "telescope_r_m": np.array([SPHERE_RADIUS_M]),
}


@pytest.fixture
def obdeect_stream():
    return pytest.importorskip("obdeect.stream")


@pytest.fixture(scope="module")
def model():
    records, files = segmented()
    return OpticalModel.from_records(records, files.__getitem__, configuration_seed=1)


def _ground_point_and_direction(count, radius_m, direction, seed):
    """Ground points (z = 0) of a beam and its float32 directions (cx, cy, cz < 0)."""
    rng = np.random.default_rng(seed)
    r = radius_m * np.sqrt(rng.random(count))
    phi = 2.0 * np.pi * rng.random(count)
    x = (r * np.cos(phi)).astype(np.float32).astype(np.float64)
    y = (r * np.sin(phi)).astype(np.float32).astype(np.float64)
    cx = np.full(count, direction[0], np.float32).astype(np.float64)
    cy = np.full(count, direction[1], np.float32).astype(np.float64)
    cz = -np.sqrt(1.0 - cx * cx - cy * cy)
    return x, y, cx, cy, cz


def plane_and_sphere_photons(model, count, direction=(0.01, -0.005), seed=4):
    """Downward photons as plane rows, and as sphere rows: the entry points of the same lines.

    The plane row is the ground point (x, y, 0) of the ray; the sphere row is where the ray enters
    the sphere of radius SPHERE_RADIUS_M centred on the origin, with its time there.
    """
    x, y, cx, cy, cz = _ground_point_and_direction(count, 0.5, direction, seed)
    ones, wavelength = np.ones(count, np.float32), np.full(count, 400.0, np.float32)
    zeros = np.zeros(count, np.float32)
    plane = {
        "x": x.astype(np.float32),
        "y": y.astype(np.float32),
        "cx": cx.astype(np.float32),
        "cy": cy.astype(np.float32),
        "t": zeros,
        "w": ones,
        "wavelength": wavelength,
    }
    b = x * cx + y * cy
    c = x * x + y * y - SPHERE_RADIUS_M**2
    s = -b - np.sqrt(b * b - c)  # Smaller root: the entry point (s < 0 for a ray from above).
    sphere = {
        "x": (x + s * cx).astype(np.float32),
        "y": (y + s * cy).astype(np.float32),
        "z": (s * cz).astype(np.float32),
        "cx": cx.astype(np.float32),
        "cy": cy.astype(np.float32),
        "cz": cz.astype(np.float32),
        "t": (s / SPEED_OF_LIGHT_M_PER_NS).astype(np.float32),
        "zem": zeros,
        "w": ones,
        "wavelength": wavelength,
    }
    return plane, sphere


def _trace(model, tmp_path, obdeect_stream, photons, run, name, **options):
    source = write_stream(tmp_path / f"{name}.in", {1: photons}, run=run)
    output = tmp_path / f"{name}.out"
    counts = trace_stream(model, source, output, seed=SEED, pointing=(0.0, 0.0), **options)
    return counts, _camera_or_empty(output, obdeect_stream), output


def _camera_or_empty(path, obdeect_stream):
    """Camera-frame columns of a stream; an empty dict when no photon was written."""
    try:
        return camera_photons(path, obdeect_stream)
    except IndexError:
        return {}


def write_sphere_chunks(path, photons, chunk, obdeect_stream):
    """A sphere-run stream of one telescope whose photon chunks hold `chunk` rows each."""
    rows = len(photons["w"])
    with obdeect_stream.Writer(str(path), "pytest") as writer:
        writer.write("run", SPHERE_RUN)
        writer.write("event", {"event_number": 1, "zenith_rad": 0.0, "azimuth_rad": 0.0})
        writer.write("telescope", {"telescope_id": 1, "frame": 1})
        for start in range(0, rows, chunk):
            part = {k: v[start : start + chunk] for k, v in photons.items()}
            writer.write_photons(1, 1, start, part)
        writer.write("telescope_end", {"telescope_id": 1, "photons": float(rows)})
        writer.write("event_end", {})
        writer.write("run_end", {})
    return path


def test_sphere_rows_trace_like_the_plane_rows_of_the_same_photons(model, tmp_path, obdeect_stream):
    plane, sphere = plane_and_sphere_photons(model, 20000)
    counts_plane, plane_out, _ = _trace(model, tmp_path, obdeect_stream, plane, RUN, "p")
    counts_sphere, sphere_out, _ = _trace(model, tmp_path, obdeect_stream, sphere, SPHERE_RUN, "s")
    assert counts_plane["photons_in"] == counts_sphere["photons_in"] == 20000
    assert counts_plane["photons_out"] > 5000
    assert counts_sphere["photons_out"] == counts_plane["photons_out"]
    # Maximum differences (float32 rounding of the start point): positions [m], directions,
    # times [ns]; the weights are identical.
    differences = {
        name: float(np.max(np.abs(plane_out[name].astype(np.float64) - sphere_out[name])))
        for name in ("x", "y", "cx", "cy", "cz", "t", "w")
    }
    assert differences["x"] < 1e-4 and differences["y"] < 1e-4
    assert differences["cx"] < 1e-6 and differences["cy"] < 1e-6 and differences["cz"] < 1e-6
    assert differences["t"] < 1e-3
    assert differences["w"] == 0.0


def test_sphere_rows_native_and_python_write_the_same_stream(model, tmp_path, obdeect_stream):
    native = _native_program()
    if native is None:
        pytest.skip("obdeect-ray-tracing-trace-stream is not installed")
    plane, sphere = plane_and_sphere_photons(model, 3000)
    horizontal = {
        **sphere,
        "cz": np.zeros(3000, np.float32),
        "cx": np.ones(3000, np.float32),
        "cy": np.zeros(3000, np.float32),
    }
    horizontal["z"] = np.full(3000, 3.0, np.float32)
    horizontal["x"] = np.full(3000, -2.0, np.float32)
    horizontal["y"] = np.linspace(-1.0, 1.0, 3000).astype(np.float32)
    model_file = tmp_path / "segmented.json"
    model.save(model_file)
    with_path = {**sphere, "path": np.full(3000, 3.0, np.float32)}
    # A plane row's path must cover its move to the input plane (about 20 m for these rows).
    plane_with_path = {**plane, "path": np.full(3000, 50.0, np.float32)}
    cases = (
        ("down", with_path, SPHERE_RUN),
        ("horizontal", {**horizontal, "path": np.full(3000, 3.0, np.float32)}, SPHERE_RUN),
        ("plane", plane_with_path, RUN),
    )
    for name, photons, run in cases:
        source = write_stream(tmp_path / f"{name}.in", {1: photons}, run=run)
        python_out = tmp_path / f"{name}.python.pst"
        native_out = tmp_path / f"{name}.native.pst"
        trace_stream(model, source, python_out, seed=SEED, pointing=(0.0, 0.0))
        argv = _native_argv(
            native,
            [str(model_file)],
            source,
            native_out,
            telescopes=None,
            pointing=(0.0, 0.0),
            transmission="full",
            threads=1,
        )
        subprocess.run(argv, check=True, capture_output=True)
        _assert_same_stream(python_out, native_out, obdeect_stream)


def test_path_is_carried_to_the_camera_frame(model, tmp_path, obdeect_stream):
    # The output path is the input path (plus the move to a plane) and the traced path: each
    # photon's path grows by at least its straight distance from the start to the focal surface.
    _, sphere = plane_and_sphere_photons(model, 2000)
    path_in = np.full(2000, 5.0, np.float32)
    counts, out, _ = _trace(
        model, tmp_path, obdeect_stream, {**sphere, "path": path_in}, SPHERE_RUN, "path"
    )
    assert counts["photons_out"] > 0
    focal = model.document["optics"]["focal_surface"]["plane_z_m"]
    # From the entry point (|r| = R) to the focal surface: at least R minus the focal distance
    # from the centre (the detected points lie within about 1 m of the axis here).
    grown = out["path"].astype(np.float64) - 5.0
    assert np.all(grown >= SPHERE_RADIUS_M - abs(focal) - 1.5)


def test_horizontal_sphere_rows_on_a_tilted_telescope_are_traced(model, tmp_path, obdeect_stream):
    # Horizontal ground photons (cz = 0) whose traced direction points down (z < 0): they start at
    # their sphere entry points and pass near the telescope axis, so they take the sequential path.
    zenith, azimuth = 30.0, 40.0
    rotation = ground_to_telescope_rotation(zenith, azimuth)
    direction = np.array([-rotation[2, 0], -rotation[2, 1], 0.0])
    direction /= np.linalg.norm(direction)
    assert (rotation @ direction)[2] < -0.4
    count = 4000
    rng = np.random.default_rng(9)
    offset = rng.uniform(-0.5, 0.5, count)  # Perpendicular to the photons' direction [m].
    perpendicular = np.array([-direction[1], direction[0], 0.0])
    # A line at perpendicular distance `offset` from the axis enters the sphere at distance
    # sqrt(R^2 - offset^2) before the axis, against d: the row (horizontal, cz = 0).
    along = np.sqrt(SPHERE_RADIUS_M**2 - offset**2)
    entry = -along[:, None] * direction + offset[:, None] * perpendicular
    photons = {
        "x": entry[:, 0].astype(np.float32),
        "y": entry[:, 1].astype(np.float32),
        "z": entry[:, 2].astype(np.float32),
        "cx": np.full(count, direction[0], np.float32),
        "cy": np.full(count, direction[1], np.float32),
        "cz": np.zeros(count, np.float32),
        "t": np.zeros(count, np.float32),
        "zem": np.zeros(count, np.float32),
        "w": np.ones(count, np.float32),
        "wavelength": np.full(count, 400.0, np.float32),
    }
    source = write_stream(tmp_path / "tilt.in", {1: photons}, run=SPHERE_RUN)
    counts = trace_stream(
        model, source, tmp_path / "tilt.out", seed=SEED, pointing=(zenith, azimuth)
    )
    assert counts["photons_in"] == count
    # The dish front faces the sky: photons from the side reach it and some are detected.
    assert counts["photons_out"] > 0


def test_start_inside_the_sphere_moving_towards_the_primary_or_the_camera(
    model, tmp_path, obdeect_stream
):
    # Rows at z = 4 m (above the dish, below the focal surface), inside the sphere.
    count = 3000
    rng = np.random.default_rng(3)
    radius = 0.4 * np.sqrt(rng.random(count))
    angle = 2.0 * np.pi * rng.random(count)
    base = {
        "x": (radius * np.cos(angle)).astype(np.float32),
        "y": (radius * np.sin(angle)).astype(np.float32),
        "z": np.full(count, 4.0, np.float32),
        "cx": np.zeros(count, np.float32),
        "cy": np.zeros(count, np.float32),
        "t": np.zeros(count, np.float32),
        "zem": np.zeros(count, np.float32),
        "w": np.ones(count, np.float32),
        "wavelength": np.full(count, 400.0, np.float32),
    }
    towards_primary = {**base, "cz": np.full(count, -1.0, np.float32)}
    towards_camera = {**base, "cz": np.full(count, 1.0, np.float32)}
    down_counts, _, _ = _trace(model, tmp_path, obdeect_stream, towards_primary, SPHERE_RUN, "down")
    up_counts, _, _ = _trace(model, tmp_path, obdeect_stream, towards_camera, SPHERE_RUN, "up")
    # Towards the primary: sequential, the dish reflects them up to the focal surface.
    assert down_counts["photons_out"] > 0.5 * count
    # Towards the camera: the first-hit path; they reach the focal surface from below.
    assert up_counts["photons_out"] == 0


def test_sphere_rows_do_not_depend_on_chunking_or_threads(model, tmp_path, obdeect_stream):
    plane, sphere = plane_and_sphere_photons(model, 5000)
    single = write_stream(tmp_path / "one.in", {1: sphere}, run=SPHERE_RUN)
    chunked = write_sphere_chunks(tmp_path / "chunked.in", sphere, 700, obdeect_stream)
    outputs = []
    for name, source, threads in (("a", single, 1), ("b", chunked, 3)):
        out = tmp_path / f"{name}.out"
        trace_stream(model, source, out, seed=SEED, pointing=(0.0, 0.0), threads=threads)
        outputs.append(camera_photons(out, obdeect_stream))
    assert outputs[0]["w"].size > 0
    for name in outputs[0]:
        assert outputs[0][name].tobytes() == outputs[1][name].tobytes(), name


def test_enclosure_rejects_a_sphere_that_does_not_enclose_the_optics(
    model, tmp_path, obdeect_stream
):
    plane, sphere = plane_and_sphere_photons(model, 10)
    small = {**RUN, "ground_representation": "sphere", "telescope_r_m": np.array([0.5])}
    with pytest.raises(ModelError, match="must enclose the optical model"):
        _trace(model, tmp_path, obdeect_stream, sphere, small, "small")


def test_enclosure_distances_name_every_element(model):
    distances = dict(enclosure_distances(model))
    assert "primary" in distances and "focal surface (plane)" in distances
    assert all(value >= 0.0 for value in distances.values())
    focal = model.document["optics"]["focal_surface"]["plane_z_m"]
    assert distances["focal surface (plane)"] == pytest.approx(abs(focal), abs=1e-9)
    assert max(distances.values()) < SPHERE_RADIUS_M


def test_sphere_rows_need_their_columns_and_a_radius(model, tmp_path, obdeect_stream):
    plane, sphere = plane_and_sphere_photons(model, 10)
    without_z = {k: v for k, v in sphere.items() if k != "z"}
    with pytest.raises(ModelError, match="lack a required column"):
        _trace(model, tmp_path, obdeect_stream, without_z, SPHERE_RUN, "noz")
    outside = {**sphere, "x": sphere["x"] + np.float32(3.0 * SPHERE_RADIUS_M)}
    with pytest.raises(ModelError, match="outside the telescope sphere"):
        _trace(model, tmp_path, obdeect_stream, outside, SPHERE_RUN, "outside")
    no_radius = {k: v for k, v in SPHERE_RUN.items() if k != "telescope_r_m"}
    with pytest.raises(ModelError, match="telescope_r_m"):
        _trace(model, tmp_path, obdeect_stream, sphere, no_radius, "norad")
    # The writer refuses unknown representations; the reader and tracer share the same rule.
    with pytest.raises(RuntimeError, match="ground_representation"):
        _trace(
            model,
            tmp_path,
            obdeect_stream,
            sphere,
            {**SPHERE_RUN, "ground_representation": "cube"},
            "badrep",
        )
    negative = {**sphere, "path": np.full(10, -1.0, np.float32)}
    with pytest.raises(Exception):  # The reader rejects a negative path (format rule).
        _trace(model, tmp_path, obdeect_stream, negative, SPHERE_RUN, "negpath")


def test_plane_runs_still_reject_horizontal_photons(model, tmp_path, obdeect_stream):
    plane, _ = plane_and_sphere_photons(model, 10)
    plane = {
        **plane,
        "cx": np.ones(10, np.float32),
        "cy": np.zeros(10, np.float32),
        "cz": np.zeros(10, np.float32),
    }
    with pytest.raises(ModelError, match="horizontal photons"):
        _trace(model, tmp_path, obdeect_stream, plane, RUN, "horizontal_plane")


def test_model_hash_is_unchanged_by_the_sphere_helpers(model):
    assert model_hash(model.document) == model.document["sha256"]
    assert minimum_input_plane_height_m(model, 0.0, 0.0) > 0.0


def test_native_program_rejects_a_sphere_that_does_not_enclose_the_optics(
    model, tmp_path, obdeect_stream
):
    native = _native_program()
    if native is None:
        pytest.skip("obdeect-ray-tracing-trace-stream is not installed")
    _, sphere = plane_and_sphere_photons(model, 10)
    small = {**RUN, "ground_representation": "sphere", "telescope_r_m": np.array([0.5])}
    source = write_stream(tmp_path / "small.in", {1: sphere}, run=small)
    model_file = tmp_path / "segmented.json"
    model.save(model_file)
    argv = _native_argv(
        native,
        [str(model_file)],
        source,
        tmp_path / "small.out",
        telescopes=None,
        pointing=(0.0, 0.0),
        transmission="full",
        threads=1,
    )
    result = subprocess.run(argv, capture_output=True, text=True)
    assert result.returncode != 0
    assert "must enclose the optical model" in result.stderr
