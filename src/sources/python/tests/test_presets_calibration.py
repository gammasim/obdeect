"""Calibration presets (flasher, led, point-source, beam-scan): run descriptions from arguments."""

import argparse
import math

import pytest
from obdeect.sources.presets import beam_scan, flasher, led, point_source

RUN = ["--observation-level", "2150", "--spectrum", "fixed:wavelength_nm=400"]


def parse(module, argv):
    parser = argparse.ArgumentParser()
    module.add_arguments(parser)
    return parser.parse_args(RUN + argv)


def drop(argv, option, count=1):
    """argv without ``option`` and its ``count`` values."""
    index = argv.index(option)
    return argv[:index] + argv[index + 1 + count :]


def assert_normalizes(description):
    """normalize() accepts the description when the config module exists."""
    config = pytest.importorskip("obdeect.sources.config")
    normalized = config.normalize(description)
    assert normalized["events"]


def close(a, b):
    return all(math.isclose(x, y, abs_tol=1e-9) for x, y in zip(a, b, strict=True))


# flasher ----------------------------------------------------------------------------------


DISH = [
    "--type",
    "dish",
    "--photons",
    "1e6",
    "2e6",
    "--events",
    "2",
    "--position",
    "0.1,0.2,1.5",
    "--camera-radius",
    "1.2",
]
SECONDARY = [
    "--type",
    "secondary",
    "--photons",
    "1e6",
    "--flasher-offset",
    "0.2",
    "--flasher-depth",
    "0.03",
    "--flasher-inclination",
    "0",
    "--mirror-camera-distance",
    "3.5",
    "--fiducial-radius",
    "1.1",
]


def test_flasher_dish_structure():
    d = flasher.build(parse(flasher, DISH))
    assert d["run"]["frame"] == "camera"
    assert d["telescopes"] == [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 1.2}]
    assert len(d["events"]) == 2  # one per photon value
    for event, photons in zip(d["events"], [1e6, 2e6], strict=True):
        assert event["repeat"] == 2
        (em,) = event["emissions"]
        assert em["kind"] == "point"
        assert em["photons"] == photons
        assert em["position_m"] == [0.1, 0.2, 1.5]
        assert em["axis"] == [0.0, 0.0, -1.0]
        assert em["to_plane_z0"] is True
        assert em["target"] == "auto"
        assert em["time_offset_ns"] == 0.0
    assert_normalizes(d)


def test_flasher_dish_rejects_missing_geometry():
    with pytest.raises(ValueError, match="--position"):
        flasher.build(parse(flasher, drop(DISH, "--position")))
    with pytest.raises(ValueError, match="--camera-radius"):
        flasher.build(parse(flasher, drop(DISH, "--camera-radius")))


def test_flasher_dish_rejects_secondary_options():
    with pytest.raises(ValueError, match="--flasher-offset applies only to --type secondary"):
        flasher.build(parse(flasher, DISH + ["--flasher-offset", "0.2"]))


def test_flasher_secondary_structure_and_positions():
    d = flasher.build(parse(flasher, SECONDARY))
    assert d["run"]["frame"] == "telescope"
    assert d["telescopes"][0]["radius_m"] == 1.1
    (event,) = d["events"]
    assert event["repeat"] == 1
    assert len(event["emissions"]) == 4
    expected = [(0.2, 0.2), (0.2, -0.2), (-0.2, 0.2), (-0.2, -0.2)]
    for em, (x, y) in zip(event["emissions"], expected, strict=True):
        assert close(em["position_m"], [x, y, 3.53])
        assert em["target"] == "auto"
        assert em["to_plane_z0"] is False
    # inclination 0 deg: altitude -90 deg, every flasher points down
    for em in event["emissions"]:
        assert close(em["axis"], [0.0, 0.0, -1.0])
    assert_normalizes(d)


def test_flasher_secondary_delays_and_fire_modes():
    d = flasher.build(parse(flasher, SECONDARY + ["--delays", "0", "10", "20", "30"]))
    assert [em["time_offset_ns"] for em in d["events"][0]["emissions"]] == [0, 10, 20, 30]

    d = flasher.build(
        parse(
            flasher, SECONDARY + ["--fire", "3", "--events", "3", "--delays", "0", "10", "20", "30"]
        )
    )
    (event,) = d["events"]
    assert event["repeat"] == 3
    (em,) = event["emissions"]
    assert em["time_offset_ns"] == 20.0
    assert em["position_m"][:2] == pytest.approx([-0.2, 0.2])

    d = flasher.build(parse(flasher, SECONDARY + ["--fire", "alternate", "--events", "5"]))
    assert len(d["events"]) == 5
    fired = [event["metadata"]["flashers"][0] for event in d["events"]]
    assert fired == [1.0, 2.0, 3.0, 4.0, 1.0]
    assert all(len(event["emissions"]) == 1 and event["repeat"] == 1 for event in d["events"])


def test_flasher_secondary_rejects_missing_geometry():
    for option in [
        "--flasher-offset",
        "--flasher-depth",
        "--flasher-inclination",
        "--mirror-camera-distance",
        "--fiducial-radius",
    ]:
        with pytest.raises(ValueError, match=option):
            flasher.build(parse(flasher, drop(SECONDARY, option)))


def test_flasher_secondary_rejects_dish_options():
    with pytest.raises(ValueError, match="--position applies only to --type dish"):
        flasher.build(parse(flasher, SECONDARY + ["--position", "0,0,1"]))


def test_flasher_requires_photons_and_observation_level():
    with pytest.raises(SystemExit):
        parse(flasher, drop(DISH, "--photons", 2))
    parser = argparse.ArgumentParser()
    flasher.add_arguments(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["--photons", "1", "--spectrum", "fixed:wavelength_nm=400"])


# led --------------------------------------------------------------------------------------


CAMERA = """\
# test camera
Rotate 0
Pixtype 1
Pixel 0 1 0.0 0.0
Pixel 1 1 10.0 0.0
Pixel 2 1 0.0 20.0
Pixel 3 1 10.0 20.0
Pixel 4 1 -10.0 0.0
Trigger of 0 [1] 2
Trigger of 2 3 [4]
"""

LED = ["--photons", "1e3", "--camera-radius", "1.2"]


@pytest.fixture
def camera_file(tmp_path):
    path = tmp_path / "camera.dat"
    path.write_text(CAMERA)
    return str(path)


def test_led_single_list_all_pixels(camera_file):
    d = led.build(parse(led, ["--camera-file", camera_file, "--use", "1 2", "--events", "2"] + LED))
    assert d["run"]["frame"] == "camera"
    (event,) = d["events"]
    assert event["repeat"] == 2
    assert event["metadata"]["pixels"] == [1.0, 2.0]
    assert [em["position_m"] for em in event["emissions"]] == [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]
    assert all(em["kind"] == "point" and em["target"] == "auto" for em in event["emissions"])
    assert all(em["angular"] == {"kind": "cone", "angle_deg": 0.0} for em in event["emissions"])
    assert_normalizes(d)


def test_led_trigger_groups(camera_file):
    d = led.build(parse(led, ["--camera-file", camera_file, "--use-trg-all"] + LED))
    assert [e["metadata"]["pixels"] for e in d["events"]] == [[0.0, 2.0], [2.0, 3.0]]
    assert d["events"][0]["emissions"][0]["position_m"] == [0.0, 0.0, 0.0]


def test_led_trigger_random_groups_are_deterministic(camera_file):
    argv = ["--camera-file", camera_file, "--use-trg-rand", "5", "--seed", "7"] + LED
    first = led.build(parse(led, argv))
    second = led.build(parse(led, argv))
    assert first == second
    assert len(first["events"]) == 5


def test_led_require_combinations(camera_file):
    d = led.build(
        parse(led, ["--camera-file", camera_file, "--use", "+0 1 2", "--require", "2"] + LED)
    )
    # required pixel 0 plus each one of 1 and 2, in list order
    assert [e["metadata"]["pixels"] for e in d["events"]] == [[0.0, 1.0], [0.0, 2.0]]


def test_led_require_random_pixels(camera_file):
    argv = [
        "--camera-file",
        camera_file,
        "--use",
        "0 1 2 3",
        "--require",
        "2",
        "--random",
        "--events",
        "4",
    ] + LED
    d = led.build(parse(led, argv))
    assert len(d["events"]) == 4
    for event in d["events"]:
        pixels = event["metadata"]["pixels"]
        assert len(pixels) == 2 and len(set(pixels)) == 2


def test_led_slaves_add_emissions(camera_file):
    d = led.build(parse(led, ["--camera-file", camera_file, "--use", "0 [4]", "--slaves"] + LED))
    (event,) = d["events"]
    assert event["metadata"]["pixels"] == [0.0]
    assert [em["position_m"][0] for em in event["emissions"]] == [0.0, -0.1]


def test_led_rejects_bad_input(camera_file):
    with pytest.raises(ValueError, match="not defined in the camera file"):
        led.build(parse(led, ["--camera-file", camera_file, "--use", "9"] + LED))
    with pytest.raises(ValueError, match="random needs"):
        led.build(parse(led, ["--camera-file", camera_file, "--use", "0 1", "--random"] + LED))
    with pytest.raises(SystemExit):
        parse(led, ["--camera-file", camera_file, "--use", "0", "--use-trg-all"] + LED)
    with pytest.raises(SystemExit):
        parse(led, ["--camera-file", camera_file, "--use", "0"] + ["--photons", "1"])


# point-source -----------------------------------------------------------------------------


POINT = [
    "--position",
    "100,0,1000",
    "--photons",
    "1e6",
    "2e6",
    "--events",
    "3",
    "--telescope-radius",
    "12.5",
    "--direction",
    "0,0,-2",
]


def test_point_source_structure():
    d = point_source.build(parse(point_source, POINT))
    assert d["run"]["frame"] == "sky"
    assert d["telescopes"] == [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 12.5}]
    assert [e["emissions"][0]["photons"] for e in d["events"]] == [1e6, 2e6]
    assert all(e["repeat"] == 3 for e in d["events"])
    em = d["events"][0]["emissions"][0]
    assert em["kind"] == "point"
    assert em["position_m"] == [100.0, 0.0, 1000.0]
    assert em["axis"] == [0.0, 0.0, -1.0]
    assert em["target"] == "auto"
    assert_normalizes(d)


def test_point_source_default_direction_and_profile():
    argv = [a for a in POINT if a not in ("--direction", "0,0,-2")]
    em = point_source.build(parse(point_source, argv))["events"][0]["emissions"][0]
    assert em["axis"] == [0.0, 0.0, -1.0]
    assert em["angular"] == {"kind": "isotropic"}


def test_point_source_rejects_zero_direction_and_missing_radius():
    with pytest.raises(ValueError, match="direction"):
        point_source.build(
            parse(point_source, drop(POINT, "--direction") + ["--direction", "0,0,0"])
        )
    with pytest.raises(SystemExit):
        parse(point_source, drop(POINT, "--telescope-radius"))


# beam-scan --------------------------------------------------------------------------------


SCAN = [
    "--camera-radius",
    "1.23",
    "--scan-start",
    "0.1",
    "0.05",
    "--scan-extent-x",
    "0.35",
    "--scan-extent-y",
    "0.25",
    "--scan-step",
    "0.1",
    "--height",
    "0.06",
    "--photons-per-position",
    "100",
]


def test_beam_scan_structure():
    d = beam_scan.build(parse(beam_scan, SCAN))
    assert d["run"]["frame"] == "sky"
    assert d["telescopes"][0]["radius_m"] == 1.23
    (event,) = d["events"]
    # x line: 0.1, 0.2, 0.3, 0.4 (< 0.45); y line from the start: 0.15, 0.25 (start once)
    positions = [em["position_m"] for em in event["emissions"]]
    assert len(positions) == 6
    assert positions[0] == pytest.approx([0.1, 0.05, 0.06])
    assert positions[4] == pytest.approx([0.1, 0.15, 0.06])
    assert positions[-1] == pytest.approx([0.1, 0.25, 0.06])
    assert all(em["target"] == "none" and em["photons"] == 100.0 for em in event["emissions"])
    assert all(em["axis"] == [0.0, 0.0, -1.0] for em in event["emissions"])
    assert_normalizes(d)


def test_beam_scan_beam_direction():
    up = beam_scan.build(parse(beam_scan, SCAN + ["--beam-zenith", "0"]))
    assert all(em["axis"] == [0.0, 0.0, 1.0] for em in up["events"][0]["emissions"])
    tilted = beam_scan.build(
        parse(beam_scan, SCAN + ["--beam-zenith", "120", "--beam-azimuth", "90"])
    )
    axis = tilted["events"][0]["emissions"][0]["axis"]
    assert axis[2] == pytest.approx(-0.5)
    assert axis[1] == pytest.approx(-math.sqrt(3) / 2)  # azimuth 90 deg is -y (east)
    assert_normalizes(up)
    with pytest.raises(ValueError, match="horizontal beams are not supported yet"):
        beam_scan.build(parse(beam_scan, SCAN + ["--beam-zenith", "90"]))


def test_beam_scan_rejects_bad_scan_and_missing_options():
    with pytest.raises(ValueError, match="step"):
        beam_scan.build(parse(beam_scan, drop(SCAN, "--scan-step") + ["--scan-step", "0"]))
    for option, count in [
        ("--scan-step", 1),
        ("--height", 1),
        ("--scan-start", 2),
        ("--camera-radius", 1),
    ]:
        with pytest.raises(SystemExit):
            parse(beam_scan, drop(SCAN, option, count))


def test_point_source_telescope_list():
    argv = [a for a in POINT if a not in ("--telescope-radius", "12.5")]
    argv += ["--telescope", "0,0,0,12", "--telescope", "100,0,2,6"]
    d = point_source.build(parse(point_source, argv))
    assert [t["position_m"] for t in d["telescopes"]] == [[0.0, 0.0, 0.0], [100.0, 0.0, 2.0]]
    assert [t["radius_m"] for t in d["telescopes"]] == [12.0, 6.0]
    assert_normalizes(d)
