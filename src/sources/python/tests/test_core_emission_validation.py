"""The native core validates normalized emissions itself (without the Python config layer).

Errors name the key of the emission, so the Python layer does not rewrite them."""

import pytest
from obdeect.sources import _core

SPECTRUM = {"kind": "fixed", "wavelength_nm": 400.0}
PULSE = {"kind": "fixed", "time_ns": 0.0}
NO_TARGET = {"kind": "none"}


def emission(**changes):
    spec = {
        "kind": "point",
        "photons": 100.0,
        "bunch_size": 1.0,
        "time_offset_ns": 0.0,
        "to_plane_z0": False,
        "spectrum": SPECTRUM,
        "pulse": PULSE,
        "target": NO_TARGET,
        "position_m": [0.0, 0.0, 10.0],
        "axis": [0.0, 0.0, -1.0],
        "angular": {"kind": "isotropic"},
    }
    spec.update(changes)
    return spec


def emit(spec, rows=10):
    return _core.emit_chunk(spec, 1, 0, 0, 0, rows)


def test_valid_emission_emits():
    chunk = emit(emission())
    assert len(chunk["w"]) == 10


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"kind": "cone"}, "unknown emission kind"),
        ({"spectrum": {"kind": "nonsense"}}, "unknown spectrum kind"),
        ({"pulse": {"kind": "nonsense"}}, "unknown pulse kind"),
        ({"angular": {"kind": "nonsense"}}, "unknown angular kind"),
        ({"target": {"kind": "nonsense"}}, "unknown target kind"),
        ({"angular": {"kind": "cone"}}, "missing"),
        ({"photons": True}, "expected a number"),
        ({"to_plane_z0": 1}, "expected true or false"),
        ({"position_m": [0.0, 0.0]}, "three components"),
        ({"axis": [0.0, 0.0, 0.0]}, "non-zero vector"),
        ({"target": {"kind": "auto"}}, "'auto' needs the telescopes"),
        ({"bunch_size": 0.5}, "bunch_size"),
        ({"time_offset_ns": float("inf")}, "finite"),
    ],
)
def test_bad_emissions_name_the_key(change, match):
    with pytest.raises(ValueError, match=match):
        emit(emission(**change))


def test_missing_key_is_named():
    spec = emission()
    del spec["spectrum"]
    with pytest.raises(ValueError, match="emission.spectrum: missing"):
        emit(spec)


def test_angular_kind_without_parameters():
    assert len(emit(emission(angular={"kind": "lambertian"}))["w"]) == 10


def test_source_is_validated_on_creation():
    with pytest.raises(ValueError, match="photons"):
        emit(emission(photons=-1.0))


def test_plane_wave_rejects_a_target():
    spec = {
        "kind": "plane_wave",
        "photons": 10.0,
        "bunch_size": 1.0,
        "time_offset_ns": 0.0,
        "to_plane_z0": False,
        "spectrum": SPECTRUM,
        "pulse": PULSE,
        "target": {"kind": "sphere", "centre_m": [0.0, 0.0, 0.0], "radius_m": 1.0},
        "towards_source": [0.0, 0.0, 1.0],
        "disc_radius_m": 1.0,
        "disc_z_m": 0.0,
        "start_z_m": 5.0,
    }
    with pytest.raises(ValueError, match="no target"):
        emit(spec)
