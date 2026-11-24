"""Inputs are used or rejected, never dropped: the explicit rules of the importer, the bounds of the
star-light disk, the convex hull of structure elements, the PSF summary and the CLI errors."""

import itertools
import math
import random

import numpy as np
import pytest
from obdeect.ray_tracing import OpticalModel
from obdeect.ray_tracing.bounds import circumradius, obscurer_top_z_m, optics_top_z_m
from obdeect.ray_tracing.cli import main
from obdeect.ray_tracing.errors import ModelError
from obdeect.ray_tracing.importers.sim_telarray import (
    TableImportError,
    parse_rpol_table,
    read_mirror_list,
    read_response,
    read_segmentation,
)
from obdeect.ray_tracing.importers.simulation_models import (
    _in_range,
    compile_model,
    convex_hull,
)
from obdeect.ray_tracing.model import DETECTED
from obdeect.ray_tracing.psf import containment, summary
from telescopes import MIRROR_LIST, dual, record, segmented


def build(records_files, **options):
    records, files = records_files
    options.setdefault("configuration_seed", 1)
    options.setdefault("propagation_group_index", 1.0)
    return compile_model(records, files.__getitem__, **options)


def test_transmission_method_zero_has_no_amplitude():
    zero = record("telescope_transmission", [0.9, 0, 0.1, 5.0, 2.0, 0.0])
    with pytest.raises(ModelError, match="method 0"):
        build(segmented(telescope_transmission=zero))
    no_amplitude = record("telescope_transmission", [0.9, 0, 0.0, 5.0, 2.0, 0.0])
    assert build(segmented(telescope_transmission=no_amplitude))["optics"]["telescope_transmission"]


def test_dual_mirror_outer_power_is_used_or_rejected():
    outer = record("telescope_transmission", [0.9, 1, 0.1, 5.0, 2.0, 0.5])
    with pytest.raises(ModelError, match="outer power"):
        build(dual(telescope_transmission=outer))
    unit = record("telescope_transmission", [0.9, 1, 0.1, 5.0, 2.0, 1.0])
    assert build(dual(telescope_transmission=unit))["optics"]["telescope_transmission"]


def test_negative_shadow_offset_is_rejected_and_zero_selects_the_automatic_position():
    diameter = record("secondary_mirror_shadow_diameter", 20.0, "cm")
    with pytest.raises(ModelError, match="must not be negative"):
        build(
            dual(
                secondary_mirror_shadow_diameter=diameter,
                secondary_mirror_shadow_offset=record("secondary_mirror_shadow_offset", -5.0, "cm"),
            )
        )
    automatic = build(
        dual(
            secondary_mirror_shadow_diameter=diameter,
            secondary_mirror_shadow_offset=record("secondary_mirror_shadow_offset", 0.0, "cm"),
        )
    )
    assert automatic["optics"]["obscurers"][0]["name"] == "secondary-shadow"


def test_parameters_without_a_no_op_default_are_required():
    with pytest.raises(ModelError, match="missing parameter camera_body_shape"):
        build(segmented(camera_body_shape=None))
    with pytest.raises(ModelError, match="missing parameter mirror_class"):
        build(segmented(mirror_class=None))
    with pytest.raises(ModelError, match="dish_shape_length must be positive"):
        build(segmented(dish_shape_length=record("dish_shape_length", 0.0, "cm")))


def test_focus_offset_with_one_zenith_term_is_rejected():
    focus = record("focus_offset", [2.0, 20.0, 3.0, 0.0], ["cm", "deg", "cm", "cm"])
    with pytest.raises(ModelError, match="zenith"):
        build(segmented(focus_offset=focus), zenith_deg=60.0)


def test_provenance_cannot_overwrite_the_compilers_keys():
    records_files = segmented()
    records, files = records_files
    with pytest.raises(ModelError, match="written by the compiler"):
        compile_model(
            records,
            files.__getitem__,
            configuration_seed=1,
            propagation_group_index=1.0,
            provenance={"zenith_deg": 3.0, "model": "x"},
        )


def test_propagation_group_index_has_no_default():
    records, files = segmented()
    with pytest.raises(TypeError):
        # Deliberately without propagation_group_index: the call must fail.
        compile_model(records, files.__getitem__, configuration_seed=1)


def test_from_records_without_group_index_warns_and_keeps_the_former_value():
    # Callers written before the argument existed (simtools) still compile, with a warning, to
    # exactly the model of an explicit 1.0.
    records, files = segmented()
    with pytest.warns(DeprecationWarning, match="propagation_group_index"):
        implicit = OpticalModel.from_records(records, files.__getitem__, configuration_seed=1)
    explicit = OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=1, propagation_group_index=1.0
    )
    assert implicit.sha256 == explicit.sha256


def test_star_start_height_is_above_every_optical_element():
    for records_files in (segmented(), dual()):
        document = build(records_files)
        optics = document["optics"]
        start = document["sampling"]["start_z_m"]
        # The bound (without its interpolation margin) and every obscurer lie 1 m below the start.
        assert start >= optics_top_z_m(optics, margin=False) + 1.0 - 1e-12
        assert all(start >= obscurer_top_z_m(o) + 1.0 - 1e-12 for o in optics["obscurers"])


def test_star_start_height_keeps_the_old_heights_as_a_floor():
    # The start is never below the vertex heights used before (a vertex may lie above the rim).
    document = build(dual())
    vertex = document["optics"]["secondary"]["asphere"]["vertex_z_m"]
    assert document["sampling"]["start_z_m"] >= vertex + 1.0


def test_reflectivity_rounding_is_clamped_and_overshoot_is_rejected():
    assert _in_range(1.0 + 1e-12) == 1.0
    assert _in_range(-1e-12) == 0.0
    assert _in_range(0.25) == 0.25
    with pytest.raises(ModelError, match="outside"):
        _in_range(1.01)


def reference_hull(vertices):
    """The face planes of the convex hull by testing every triple against every vertex (O(n^4))."""

    def sub(a, b):
        return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    scale = max(math.dist(a, b) for a in vertices for b in vertices)
    tolerance = 1e-9 * scale
    planes = []
    for i, j, k in itertools.combinations(range(len(vertices)), 3):
        a, b, c = vertices[i], vertices[j], vertices[k]
        e, f = sub(b, a), sub(c, a)
        normal = [
            e[1] * f[2] - e[2] * f[1],
            e[2] * f[0] - e[0] * f[2],
            e[0] * f[1] - e[1] * f[0],
        ]
        length = math.sqrt(dot(normal, normal))
        if length <= 1e-12 * scale * scale:
            continue
        normal = [v / length for v in normal]
        offset = dot(normal, a)
        sides = [dot(normal, v) - offset for v in vertices]
        if all(s <= tolerance for s in sides):
            plane = [*normal, offset]
        elif all(s >= -tolerance for s in sides):
            plane = [-normal[0], -normal[1], -normal[2], -offset]
        else:
            continue
        if not any(
            dot(plane[:3], q[:3]) > 1 - 1e-12 and abs(plane[3] - q[3]) <= tolerance for q in planes
        ):
            planes.append(plane)
    return planes


def test_convex_hull_is_the_reference_hull_in_the_same_order():
    cube = [[x, y, z] for x, y, z in itertools.product((0.0, 1.0), repeat=3)]
    grid = [[float(x), float(y), 0.0] for x, y in itertools.product(range(3), repeat=2)]
    grid += [[1.0, 1.0, 2.0], [0.5, 0.5, 3.0]]
    generator = random.Random(4)
    cloud = [[generator.uniform(-1, 1) for _ in range(3)] for _ in range(12)]
    for vertices in (cube, grid, cloud):
        assert convex_hull(vertices) == reference_hull(vertices)


def test_mirror_list_panel_ids_are_labels_and_other_columns_are_rejected():
    facets = read_mirror_list(MIRROR_LIST + "10 20 30 -1500 1 -40 7\n")
    assert facets[-1]["panel_id"] == 7
    assert facets[0]["panel_id"] is None
    for line, message in (
        ("10 20 30 -1500 1 -40 7.5\n", "panel id must be an integer"),
        ("10 20 30 -1500 1 -40 7 9\n", "expected 6 or 7 columns"),
        ("10 20 30 -1500 1\n", "expected 6 or 7 columns"),
    ):
        with pytest.raises(ModelError, match=message):
            read_mirror_list(MIRROR_LIST + line)


def test_ring_segments_have_consecutive_ids():
    segments = read_segmentation("ring 3 100 200 120 0 2\n")
    assert [s["id"] for s in segments] == [0, 1, 2]


def test_empty_response_file_is_a_model_error():
    with pytest.raises(ModelError, match="is empty"):
        read_response("\n  \n", "mirror_reflectivity")


def test_rpol_options_are_used_or_rejected():
    table = "#@RPOL@ 2 scheme=2\n0 45\n100 0.5 0.6\n200 0.4 0.3\n"
    with pytest.raises(TableImportError, match="bilinear"):
        parse_rpol_table(table)
    assert parse_rpol_table("#@RPOL@ 2 scheme=1\n0 45\n100 0.5 0.6\n200 0.4 0.3\n")["x"] == [
        100.0,
        200.0,
    ]
    with pytest.raises(TableImportError, match="must not be zero"):
        parse_rpol_table("#@RPOL@ 1 xscale=0\n1 0.5\n2 0.6\n")
    with pytest.raises(TableImportError, match="unsupported RPOL option"):
        parse_rpol_table("#@RPOL@ 1 zxmax=5\n1 0.5\n2 0.6\n")
    assert parse_rpol_table("#@RPOL@ 1 zxmin=0 verbose\n1 0.5\n2 0.6\n")["x"] == [1.0, 2.0]


def test_structure_corners_need_every_coordinate():
    structure = """# %ECSV 1.0
# ---
# coordinate_frame: primary_mirror_frame
# datatype:
# - {name: id, datatype: string}
# - {name: group, datatype: string}
# - {name: shape, datatype: string}
# - {name: diameter, unit: m, datatype: float64}
# - {name: x1, unit: m, datatype: float64}
# - {name: z1, unit: m, datatype: float64}
# - {name: x2, unit: m, datatype: float64}
# - {name: z2, unit: m, datatype: float64}
id group shape diameter x1 z1 x2 z2
mast structure cylinder 0.2 0 0 1 2
"""
    with pytest.raises(ModelError, match="lacks the columns y1, y2"):
        build(segmented(telescope_structure=structure))


def test_hull_planes_of_a_cube_bound_it():
    planes = convex_hull([[x, y, z] for x, y, z in itertools.product((0.0, 1.0), repeat=3)])
    assert len(planes) == 6
    assert all(abs(math.hypot(p[0], p[1], p[2]) - 1.0) < 1e-12 for p in planes)


def test_containment_needs_positive_weight_and_a_valid_fraction():
    x = np.array([0.0, 1.0, 2.0])
    with pytest.raises(ValueError, match="sum to zero"):
        containment(x, x, np.zeros(3))
    with pytest.raises(ValueError, match="fraction"):
        containment(x, x, fraction=1.5)
    assert containment(x, x)["radius_m"] >= 0.0


def test_detected_only_summary_needs_the_sample_size_and_reports_zero_weight():
    arrivals = {
        "status": np.array([DETECTED]),
        "weight": np.array([0.5]),
        "position_m": np.zeros((1, 3)),
        "index": np.array([3]),
    }
    with pytest.raises(ValueError, match="sampled"):
        summary(arrivals, 2.0)
    result = summary(arrivals, 2.0, sampled=10)
    assert result["photons"] == 10
    assert result["effective_area_m2"] == pytest.approx(0.5 / 10 * 2.0)
    zero = {k: v for k, v in arrivals.items() if k != "index"} | {"weight": np.zeros(1)}
    no_psf = summary(zero, 2.0)
    assert "radius_m" not in no_psf and "d80_m" not in no_psf
    assert no_psf["detected_weight"] == 0.0


def test_plate_outline_follows_its_orientation():
    plot = pytest.importorskip("obdeect.ray_tracing.plot")
    square = np.array([[1.0, 0.0], [0.0, 1.0]])
    points = plot._in_plane(square, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0])
    assert points[0] == pytest.approx([0.0, 1.0, 0.0])
    assert points[1] == pytest.approx([-1.0, 0.0, 0.0])


def test_circumradius_rejects_unknown_apertures():
    assert circumradius("square", 2.0) == pytest.approx(math.sqrt(2.0))
    with pytest.raises(ModelError, match="unknown aperture"):
        circumradius("triangle", 2.0)


def test_model_positions_and_pools():
    from obdeect.ray_tracing.model import _positions

    assert _positions([5, 3], [3, 5, 7]).tolist() == [1, 0]
    with pytest.raises(ModelError, match="not an obscurer"):
        _positions([4], [3, 5, 7])
    records, files = segmented()
    model = OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=1, propagation_group_index=1.0
    )
    assert model._shared_pool(2) is model._shared_pool(2)
    with pytest.raises(ModelError, match="block_size"):
        next(
            model.trace_emitter_blocks(
                1.0, seed=1, position_m=(0, 0, 0), intensity=([0, 90], [1, 1]), block_size=0
            )
        )


def test_cli_reports_invalid_json_as_an_error(tmp_path, capsys):
    path = tmp_path / "model.json"
    path.write_text("{", encoding="utf-8")
    assert main(["trace", "--model", str(path), "--seed", "1"]) == 1
    assert "error" in capsys.readouterr().err
