"""Run description: loading, normalization, defaults and error paths (pure Python)."""

import copy
import json

import pytest
from obdeect.sources import config, run

FULL_TOML = """
[run]
number = 2
seed = 7
observation_level_m = 2150.0
atmosphere_id = 3
frame = "telescope"
metadata = { note = "x", level = 1, scale = 0.5, ids = [1, 2.5] }

[[telescopes]]
id = 5
position_m = [10.0, 0.0, 0.0]
radius_m = 12.5

[[telescopes]]
radius_m = 6.0

[[events]]
repeat = 2
metadata = { tag = "a" }

[[events.emissions]]
kind = "point"
photons = 1e6
bunch_size = 4.0
time_offset_ns = 1.5
to_plane_z0 = true
position_m = [0.0, 0.0, 1.0]
axis = [0.0, 0.0, -1.0]
angular = { kind = "gauss", sigma_deg = 3.0 }
spectrum = { kind = "uniform", low_nm = 300.0, high_nm = 600.0 }
pulse = { kind = "gauss", sigma_ns = 1.0 }
target = { kind = "disc", centre_m = [0.0, 0.0, 0.0], normal = [0.0, 0.0, 1.0], radius_m = 1.2 }

[[events.emissions]]
kind = "line"
photons = 2e5
start_m = [0.0, 0.0, 5000.0]
direction = [0.0, 0.0, -1.0]
length_m = 1000.0
beta = 0.9
angular = { kind = "table", angle_deg = [0.0, 10.0], value = [1.0, 0.5] }
profile = { distance_m = [0.0, 1000.0], density = [1.0, 2.0] }
spectrum = { kind = "table", wavelength_nm = [300.0, 600.0], density = [1.0, 0.5] }
pulse = { kind = "tophat", width_ns = 2.5 }
target = "auto"

[[events.emissions]]
kind = "plane_wave"
photons = 1e4
towards_source = [0.0, 0.0, 1.0]
distance_m = 1e5
disc_radius_m = 30.0
spectrum = { kind = "inverse_square", low_nm = 200.0, high_nm = 800.0 }
pulse = { kind = "exponential", decay_ns = 3.0 }
target = "none"

[[events.emissions]]
kind = "sky_background"
photons = 0.0
radius_m = 2000.0
spectrum = { kind = "fixed", wavelength_nm = 400.0 }
pulse = { kind = "complex", flat_ns = 1.0, gauss_rms_ns = 0.5, decay_ns = 2.0 }
target = { kind = "sphere", centre_m = [0.0, 0.0, 0.0], radius_m = 50.0 }
"""


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _minimal():
    return {
        "run": {"observation_level_m": 2150.0},
        "telescopes": [{"radius_m": 12.5}],
        "events": [
            {
                "emissions": [
                    {
                        "kind": "point",
                        "photons": 1e6,
                        "position_m": [0.0, 0.0, 1.0],
                        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
                    },
                    {
                        "kind": "point",
                        "photons": 1e6,
                        "position_m": [0.0, 0.0, 1.0],
                        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
                    },
                ]
            }
        ],
    }


def test_full_toml_round_trip(tmp_path):
    raw = config.load(_write(tmp_path, "run.toml", FULL_TOML))
    normalized = config.normalize(raw)

    assert normalized["run"] == {
        "number": 2,
        "seed": 7,
        "observation_level_m": 2150.0,
        "atmosphere_id": 3,
        "frame": "telescope",
        "metadata": {"note": "x", "level": 1, "scale": 0.5, "ids": [1.0, 2.5]},
    }
    assert normalized["telescopes"] == [
        {"id": 5, "position_m": [10.0, 0.0, 0.0], "radius_m": 12.5},
        {"id": 2, "position_m": [0.0, 0.0, 0.0], "radius_m": 6.0},
    ]
    point, line, plane, sky = normalized["events"][0]["emissions"]
    assert point["kind"] == "point"
    assert point["target"] == {
        "kind": "disc",
        "centre_m": [0.0, 0.0, 0.0],
        "normal": [0.0, 0.0, 1.0],
        "radius_m": 1.2,
    }
    assert line["target"] == {"kind": "auto"}
    assert line["profile"] == {"distance_m": [0.0, 1000.0], "density": [1.0, 2.0]}
    assert plane["distance_m"] == 1e5
    assert sky["target"] == {"kind": "sphere", "centre_m": [0.0, 0.0, 0.0], "radius_m": 50.0}

    # Normalized form is a fixed point: it is valid input and gives itself back.
    assert config.normalize(normalized) == normalized


def test_json_gives_same_result(tmp_path):
    raw = config.load(_write(tmp_path, "run.toml", FULL_TOML))
    json_path = _write(tmp_path, "run.json", json.dumps(raw))
    assert config.normalize(config.load(json_path)) == config.normalize(raw)


def test_load_rejects_other_suffixes(tmp_path):
    path = _write(tmp_path, "run.yaml", "run: {}")
    with pytest.raises(ValueError, match=r"unsupported suffix '.yaml'"):
        config.load(path)


def test_defaults_are_filled():
    normalized = config.normalize(_minimal())

    assert normalized["run"] == {
        "number": 1,
        "seed": 1,
        "observation_level_m": 2150.0,
        "atmosphere_id": -1,
        "frame": "sky",
        "metadata": {},
    }
    assert normalized["telescopes"] == [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 12.5}]
    event = normalized["events"][0]
    assert event["repeat"] == 1 and event["metadata"] == {}
    emission = event["emissions"][0]
    assert emission == {
        "kind": "point",
        "photons": 1e6,
        "bunch_size": 1.0,
        "time_offset_ns": 0.0,
        "to_plane_z0": False,
        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
        "pulse": {"kind": "fixed", "time_ns": 0.0},
        "target": {"kind": "none"},
        "position_m": [0.0, 0.0, 1.0],
        "axis": [0.0, 0.0, 1.0],
        "angular": {"kind": "isotropic"},
    }


def test_line_and_sky_background_defaults():
    config_dict = _minimal()
    config_dict["events"][0]["emissions"] = [
        {
            "kind": "line",
            "photons": 1,
            "start_m": [0, 0, 1],
            "direction": [0, 0, -1],
            "length_m": 10,
            "spectrum": {"kind": "fixed", "wavelength_nm": 400},
        },
        {
            "kind": "sky_background",
            "photons": 1,
            "radius_m": 5,
            "spectrum": {"kind": "fixed", "wavelength_nm": 400},
        },
    ]
    line, sky = config.normalize(config_dict)["events"][0]["emissions"]
    assert line["beta"] == 1.0
    assert line["angular"] == {"kind": "cone", "angle_deg": 0.0}
    assert "profile" not in line
    assert sky["centre_m"] == [0.0, 0.0, 0.0]
    assert sky["axis"] == [0.0, 0.0, -1.0]
    assert sky["viewcone_deg"] == 90.0


def test_every_profile_kind_is_accepted():
    profiles = {
        "spectrum": [
            {"kind": "fixed", "wavelength_nm": 400.0},
            {"kind": "uniform", "low_nm": 300.0, "high_nm": 600.0},
            {"kind": "inverse_square", "low_nm": 300.0, "high_nm": 600.0},
            {"kind": "table", "wavelength_nm": [300.0, 600.0], "density": [1.0, 1.0]},
        ],
        "pulse": [
            {"kind": "fixed", "time_ns": 1.0},
            {"kind": "tophat", "width_ns": 2.0},
            {"kind": "gauss", "sigma_ns": 1.0},
            {"kind": "exponential", "decay_ns": 1.0},
            {"kind": "complex", "flat_ns": 1.0, "gauss_rms_ns": 1.0, "decay_ns": 1.0},
            {"kind": "table", "time_ns": [0.0, 1.0], "rate": [1.0, 2.0]},
        ],
        "angular": [
            {"kind": "isotropic"},
            {"kind": "cone", "angle_deg": 5.0},
            {"kind": "filled_cone", "angle_deg": 5.0},
            {"kind": "lambertian"},
            {"kind": "gauss", "sigma_deg": 2.0},
            {"kind": "rayleigh"},
            {"kind": "table", "angle_deg": [0.0, 90.0], "value": [1.0, 0.0]},
        ],
    }
    for section, variants in profiles.items():
        for variant in variants:
            config_dict = _minimal()
            emission = config_dict["events"][0]["emissions"][0]
            if section == "spectrum":
                emission["spectrum"] = variant
            elif section == "pulse":
                emission["pulse"] = variant
            else:
                emission["angular"] = variant
            emission = config.normalize(config_dict)["events"][0]["emissions"][0]
            assert emission[section]["kind"] == variant["kind"]


def _set(*path_and_value):
    """Mutator: set config[path...] = value (path items are keys or list indices)."""
    *path, value = path_and_value

    def mutate(cfg):
        node = cfg
        for key in path[:-1]:
            node = node[key]
        node[path[-1]] = value

    return mutate


def _del(*path):
    def mutate(cfg):
        node = cfg
        for key in path[:-1]:
            node = node[key]
        del node[path[-1]]

    return mutate


ERROR_CASES = [
    (_del("run", "observation_level_m"), "run.observation_level_m: required"),
    (_del("run"), "run: required"),
    (_del("events"), "events: required"),
    (_set("events", []), "events: needs at least one entry"),
    (_set("telescopes", []), "telescopes: needs at least one entry"),
    (_set("run", "foo", 1), "run.foo: unknown key"),
    (_set("foo", 1), "foo: unknown key"),
    (_set("run", "number", 1.0), "run.number: expected an int, got float"),
    (_set("run", "number", True), "run.number: expected an int, got bool"),
    (_set("run", "number", 0), "run.number: must be >= 1"),
    (_set("run", "seed", -1), "run.seed: must be in 0..9223372036854775807"),
    (_set("run", "seed", 2**63), "run.seed: must be in 0..9223372036854775807"),
    (_set("run", "frame", "earth"), "run.frame: must be one of sky, telescope, camera, ground"),
    (_set("run", "observation_level_m", float("nan")), "run.observation_level_m: must be finite"),
    (
        _set("run", "metadata", {"note": True}),
        "run.metadata.note: expected int, float, str or list of numbers, got bool",
    ),
    (_set("run", "metadata", {"d": [1, "x"]}), "run.metadata.d[1]: expected a number, got str"),
    (_set("telescopes", 0, "radius_m", 0.0), "telescopes[0].radius_m: must be > 0"),
    (
        _set("telescopes", 0, "position_m", [1.0, 2.0]),
        "telescopes[0].position_m: expected a list of 3 numbers",
    ),
    (_set("telescopes", 0, "id", 0), "telescopes[0].id: must be >= 1"),
    (
        lambda cfg: cfg["telescopes"].append({"id": 1, "radius_m": 5.0}),
        "telescopes[1].id: duplicate id 1 (also telescopes[0])",
    ),
    (_set("events", 0, "repeat", -1), "events[0].repeat: must be >= 0"),
    (
        _set("events", 0, "emissions", 0, "photons", -5.0),
        "events[0].emissions[0].photons: must be >= 0",
    ),
    (
        _set("events", 0, "emissions", 0, "photons", True),
        "events[0].emissions[0].photons: expected a number, got bool",
    ),
    (
        _set("events", 0, "emissions", 0, "bunch_size", 0.5),
        "events[0].emissions[0].bunch_size: must be >= 1",
    ),
    (
        _set("events", 0, "emissions", 0, "to_plane_z0", 1),
        "events[0].emissions[0].to_plane_z0: expected true or false, got int",
    ),
    (_del("events", 0, "emissions", 0, "kind"), "events[0].emissions[0].kind: required"),
    (
        _set("events", 0, "emissions", 0, "kind", "cones"),
        "events[0].emissions[0].kind: unknown kind 'cones' "
        "(one of point, line, plane_wave, sky_background)",
    ),
    (_del("events", 0, "emissions", 0, "spectrum"), "events[0].emissions[0].spectrum: required"),
    (
        _set("events", 0, "emissions", 1, "angular", {"kind": "cones"}),
        "events[0].emissions[1].angular.kind: unknown kind 'cones' (one of isotropic, cone, "
        "filled_cone, lambertian, gauss, rayleigh, table)",
    ),
    (
        _set("events", 0, "emissions", 1, "angular", {"kind": "cone"}),
        "events[0].emissions[1].angular.angle_deg: required",
    ),
    (
        _set("events", 0, "emissions", 1, "angular", {"kind": "cone", "angle_deg": 5, "x": 1}),
        "events[0].emissions[1].angular.x: unknown key",
    ),
    (
        _set(
            "events",
            0,
            "emissions",
            0,
            "spectrum",
            {"kind": "table", "wavelength_nm": [1, 2, 3], "density": [1, 2]},
        ),
        "events[0].emissions[0].spectrum: wavelength_nm and density differ in length",
    ),
    (
        _set(
            "events",
            0,
            "emissions",
            0,
            "spectrum",
            {"kind": "table", "wavelength_nm": [1], "density": [1]},
        ),
        "events[0].emissions[0].spectrum.wavelength_nm: expected a list of at least 2 numbers",
    ),
    (
        _set(
            "events",
            0,
            "emissions",
            0,
            "spectrum",
            {"kind": "table", "wavelength_nm": [1, float("inf")], "density": [1, 2]},
        ),
        "events[0].emissions[0].spectrum.wavelength_nm[1]: must be finite",
    ),
    (
        _set("events", 0, "emissions", 0, "target", "disk"),
        "events[0].emissions[0].target.kind: unknown kind 'disk' (one of none, auto, disc, sphere)",
    ),
    (
        _set(
            "events",
            0,
            "emissions",
            0,
            "target",
            {"kind": "disc", "centre_m": [0, 0, 0], "normal": [0, 0, 1]},
        ),
        "events[0].emissions[0].target.radius_m: required",
    ),
    (
        _set(
            "events",
            0,
            "emissions",
            0,
            "target",
            {"kind": "sphere", "centre_m": [0, 0, 0], "radius_m": 0},
        ),
        "events[0].emissions[0].target.radius_m: must be > 0",
    ),
    (
        _set("events", 0, "emissions", 0, "target", {"kind": "none", "radius_m": 1}),
        "events[0].emissions[0].target.radius_m: unknown key",
    ),
]


@pytest.mark.parametrize(("mutate", "message"), ERROR_CASES)
def test_errors_name_the_key_path(mutate, message):
    cfg = _minimal()
    mutate(cfg)
    with pytest.raises(ValueError) as excinfo:
        config.normalize(cfg)
    assert str(excinfo.value) == message


def test_input_is_not_mutated(tmp_path):
    raw = config.load(_write(tmp_path, "run.toml", FULL_TOML))
    before = copy.deepcopy(raw)
    normalized = config.normalize(raw)
    assert raw == before
    normalized["run"]["metadata"]["note"] = "changed"
    normalized["events"][0]["emissions"][1]["profile"]["density"][0] = 99.0
    assert raw == before


def test_ground_frame_and_exactly_one_start_key():
    description = {
        "run": {"observation_level_m": 2150.0, "frame": "ground"},
        "telescopes": [{"id": 1, "radius_m": 10.0}],
        "events": [
            {
                "emissions": [
                    {
                        "kind": "plane_wave",
                        "photons": 1.0,
                        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
                        "towards_source": [1.0, 0.0, 0.0],
                        "disc_radius_m": 2.0,
                        "start_distance_m": 30.0,
                    }
                ]
            }
        ],
    }
    normalized = config.normalize(description)
    assert normalized["run"]["frame"] == "ground"
    emission = normalized["events"][0]["emissions"][0]
    assert emission["start_distance_m"] == 30.0
    assert "start_z_m" not in emission
    # The start-key rules are the core's (checked when the run is checked, not normalized).
    both = copy.deepcopy(description)
    both["events"][0]["emissions"][0]["start_z_m"] = 30.0
    with pytest.raises(ValueError, match="not both"):
        run.check(both)
    # Neither key: the start plane is z = 0 (vertical directions); horizontal needs the distance.
    neither = copy.deepcopy(description)
    del neither["events"][0]["emissions"][0]["start_distance_m"]
    neither["events"][0]["emissions"][0]["towards_source"] = [0.0, 0.0, 1.0]
    assert "start_z_m" not in config.normalize(neither)["events"][0]["emissions"][0]
    horizontal = copy.deepcopy(neither)
    horizontal["events"][0]["emissions"][0]["towards_source"] = [1.0, 0.0, 0.0]
    with pytest.raises(ValueError, match="start_distance_m"):
        run.check(horizontal)
