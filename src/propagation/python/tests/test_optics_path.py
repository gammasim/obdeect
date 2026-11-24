"""In-telescope path length from an optical model (synthetic test geometry, not CTAO)."""

import pytest
from obdeect.propagation.optics_path import _path_m, _vertex_z_m, in_telescope_path_m
from obdeect.ray_tracing.model import OpticalModel

MIRROR_LIST = """# x y diameter focal shape z (z = 0: automatic height)
   0.0     0.0  100.0  0.0  0  0.0
 125.0     0.0  100.0  0.0  0  0.0
-125.0     0.0  100.0  0.0  0  0.0
"""


def record(name, value, unit=None, *, file=False):
    return {
        "parameter": name,
        "parameter_version": "1.0.0",
        "value": value,
        "unit": unit,
        "file": file,
    }


def dual_records(*, mirror_offset_cm=0.0):
    """Dual mirror: primary vertex 0, secondary at 3.0 m, focus at 2.5 m (before the offset)."""
    cm = ["cm"]
    return {
        "mirror_class": record("mirror_class", 2),
        "focal_length": record("focal_length", 250.0, "cm"),
        "mirror_offset": record("mirror_offset", mirror_offset_cm, "cm"),
        "primary_mirror_parameters": record("primary_mirror_parameters", [0.0], cm),
        "primary_mirror_diameter": record("primary_mirror_diameter", 400.0, "cm"),
        "primary_mirror_hole_diameter": record("primary_mirror_hole_diameter", 0.0, "cm"),
        "secondary_mirror_parameters": record("secondary_mirror_parameters", [300.0], cm),
        "secondary_mirror_diameter": record("secondary_mirror_diameter", 100.0, "cm"),
        "secondary_mirror_hole_diameter": record("secondary_mirror_hole_diameter", 0.0, "cm"),
        "focal_surface_parameters": record("focal_surface_parameters", [250.0], cm),
        "camera_body_diameter": record("camera_body_diameter", 50.0, "cm"),
        "camera_body_shape": record("camera_body_shape", 2),
        "camera_depth": record("camera_depth", 0.0, "cm"),
    }


def segmented_records():
    """Single mirror with facets: no primary vertex is stored."""
    return {
        "mirror_class": record("mirror_class", 0),
        "focal_length": record("focal_length", 800.0, "cm"),
        "mirror_offset": record("mirror_offset", 0.0, "cm"),
        "dish_shape_length": record("dish_shape_length", 800.0, "cm"),
        "parabolic_dish": record("parabolic_dish", False),
        "mirror_focal_length": record("mirror_focal_length", 0.0, "cm"),
        "mirror_list": record("mirror_list", "mirrors.dat", file=True),
        "camera_body_diameter": record("camera_body_diameter", 60.0, "cm"),
        "camera_body_shape": record("camera_body_shape", 2),
        "camera_depth": record("camera_depth", 0.0, "cm"),
    }


def build(records, files=None):
    files = files or {}
    return OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=1, propagation_group_index=1.0
    )


def test_path_single_mirror_formula():
    # Focal length from the primary vertex to the focus: |f0 - p0|.
    assert _path_m(0.5, None, 8.5) == pytest.approx(8.0)
    assert _path_m(8.5, None, 0.5) == pytest.approx(8.0)


def test_path_dual_mirror_formula():
    # |p0 - s0| + |s0 - f0| with the secondary between primary and focus.
    assert _path_m(0.0, 3.0, 2.5) == pytest.approx(3.5)
    # Secondary beyond the focus: the formula uses absolute distances.
    assert _path_m(0.0, 3.0, 4.0) == pytest.approx(4.0)


@pytest.mark.parametrize("mirror_offset_cm", [0.0, -50.0])
def test_dual_mirror_model(mirror_offset_cm):
    # Primary vertex 0 (shifted by the mirror offset), secondary 3.0 m, focus 2.5 m.
    # The mirror offset moves all three positions equally, so the path is 3.5 m.
    model = build(dual_records(mirror_offset_cm=mirror_offset_cm))
    assert in_telescope_path_m(model) == pytest.approx(3.5)


def test_dual_mirror_model_from_file(tmp_path):
    model = build(dual_records())
    path = tmp_path / "model.json"
    model.save(path)
    assert in_telescope_path_m(str(path)) == pytest.approx(in_telescope_path_m(model))
    assert in_telescope_path_m(path) == pytest.approx(3.5)


@pytest.mark.parametrize("mirror_offset_cm", [0.0, -50.0, 30.0])
def test_segmented_single_mirror_model(mirror_offset_cm):
    # Compiled single mirror with facets: the dish vertex is at -mirror_offset, the focus at the
    # focal length above it, so the path is the focal length (8 m) for any offset.
    records = segmented_records()
    records["mirror_offset"] = record("mirror_offset", mirror_offset_cm, "cm")
    model = build(records, {"mirror_list": MIRROR_LIST})
    assert model.document["optics"]["primary"]["vertex_z_m"] == pytest.approx(
        -mirror_offset_cm / 100.0
    )
    assert in_telescope_path_m(model) == pytest.approx(8.0)


def test_segmented_model_from_file(tmp_path):
    model = build(segmented_records(), {"mirror_list": MIRROR_LIST})
    path = tmp_path / "segmented.json"
    model.save(path)
    assert in_telescope_path_m(path) == pytest.approx(8.0)


def test_segmented_mirror_without_stored_vertex_is_rejected():
    # Models compiled before the vertex was stored lack the key (the loader checks the document's
    # sha256, so such a file is rejected there; this is the path's own check): not guessed.
    with pytest.raises(ValueError, match="compile the model again"):
        _vertex_z_m({"facets": [{"id": 0}]}, "primary")
    assert _vertex_z_m({"facets": [], "vertex_z_m": -0.5}, "primary") == pytest.approx(-0.5)


def test_pathlen_command_prints_metres(tmp_path, capsys):
    from obdeect.propagation.optics_path import main

    model = build(segmented_records(), {"mirror_list": MIRROR_LIST})
    path = tmp_path / "segmented.json"
    model.save(path)
    assert main([str(path)]) == 0
    assert float(capsys.readouterr().out) == pytest.approx(8.0)


def test_pathlen_command_reports_missing_file(tmp_path, capsys):
    from obdeect.propagation.optics_path import main

    assert main([str(tmp_path / "missing.json")]) == 1
    assert capsys.readouterr().err.startswith("obdeect propagation pathlen: ")


def test_model_without_optics_section_is_not_an_optical_model():
    model = build(segmented_records(), {"mirror_list": MIRROR_LIST})
    document = {key: value for key, value in model.document.items() if key != "optics"}
    model._document = document  # the compiled document, without its optics section
    with pytest.raises(ValueError, match="has no 'optics' section: it is not an optical model"):
        in_telescope_path_m(model)


def test_exported_from_the_package():
    import obdeect.propagation as propagation

    assert propagation.in_telescope_path_m is in_telescope_path_m


def test_rejects_unsupported_argument():
    with pytest.raises(TypeError):
        in_telescope_path_m(3.5)
