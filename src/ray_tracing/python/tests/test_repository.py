"""Reading full and patch productions from simulation-models checkouts."""

import json

import pytest
from obdeect.ray_tracing.errors import ModelError
from obdeect.ray_tracing.repository import read_telescope
from telescopes import segmented, write_simulation_models


def test_patch_production_overlays_its_history(tmp_path):
    root = tmp_path / "simulation-models"
    records, files = segmented()
    write_simulation_models(root, "TEST-design", records, files)
    patch = root / "productions" / "9.9.10"
    patch.mkdir()
    (patch / "info.yml").write_text(
        'model_version: "9.9.10"\nmodel_update: "patch_update"\n'
        'model_version_history:\n  - "9.9.9"\nmodel_status: "production"\n'
    )
    manifest = {
        "model_version": "9.9.10",
        "parameters": {"TEST-design": {"focal_length": "2.0.0"}},
        "deprecated_parameters": ["camera_depth"],
    }
    (patch / "TEST-design.json").write_text(json.dumps(manifest))
    changed_focal = records["focal_length"] | {"parameter_version": "2.0.0", "value": 810.0}
    directory = root / "model_parameters" / "TEST-design" / "focal_length"
    (directory / "focal_length-2.0.0.json").write_text(json.dumps(changed_focal))

    patched = read_telescope(root, "TEST-design", "9.9.10")
    assert patched.records["focal_length"]["value"] == 810.0
    assert "camera_depth" not in patched.records
    assert set(patched.records) == set(records) - {"camera_depth"}
    assert read_telescope(root, "TEST-design", "9.9.9").records["focal_length"]["value"] == 800.0
    with pytest.raises(ModelError, match="no production table MISSING in 9.9.10"):
        read_telescope(root, "MISSING", "9.9.10")
