"""Beam presets (laser-beam, lidar, night-sky) and the beam geometry they use."""

import argparse
import math

import pytest
from obdeect.sources.presets import _beam, laser_beam, lidar, night_sky

Z_AXIS = (0.0, 0.0, 1.0)
DOWN = (0.0, 0.0, -1.0)

LASER = (
    "--observation-level 2150 --spectrum fixed:wavelength_nm=355 "
    "--laser-position 1000,0,0 --laser-zenith 30 --laser-azimuth 180 --photons-per-m 1e4 "
    "--telescope-zenith 0 --telescope-azimuth 0 --telescope-radius 12.5 --fov-deg 15"
)
LIDAR = (
    "--observation-level 2150 --satellite-position 1000,0,4000 --beam-zenith 180 "
    "--beam-azimuth 0 --divergence 0.1 --pulse-energy 1e-6 --wavelength-nm 532 "
    "--telescope-zenith 0 --telescope-azimuth 0 --telescope-radius 12.5 --fov-deg 80"
)
NIGHT_SKY = (
    "--observation-level 2150 --spectrum table:{table} --photons 1e6 --viewcone-deg 10 "
    "--telescope-position 0,0,0 --telescope-radius 12.5"
)


def parse(module, text):
    """Namespace of a preset from its option text."""
    parser = argparse.ArgumentParser()
    module.add_arguments(parser)
    return parser.parse_args(text.split())


def without(text, option):
    """Option text with one ``--option VALUE`` pair removed."""
    tokens = text.split()
    index = tokens.index(option)
    return " ".join(tokens[:index] + tokens[index + 2 :])


def emission(description):
    """The single emission of the first event of a run description."""
    emissions = description["events"][0]["emissions"]
    assert len(emissions) == 1
    return emissions[0]


def on_cone_surface(point, axis_tan):
    """True if the point lies on the cone of tan(half angle) ``axis_tan`` about the z axis."""
    x, y, z = point
    return math.hypot(x, y) == pytest.approx(axis_tan * z, rel=1e-9)


# Geometry of the beam in the field-of-view cone (_beam).


def test_beam_along_cone_axis_is_inside_from_apex_to_start():
    enter, length = _beam.cone_path((0.0, 0.0, 100.0), DOWN, Z_AXIS, 10.0)
    assert enter == pytest.approx(0.0)
    assert length == pytest.approx(100.0)


def test_horizontal_beam_crosses_cone_at_expected_interval():
    # Cone of half angle 45 deg about z: |x| <= z at z = 10. The beam at x = -100 enters at
    # x = -10 (t = 90) and leaves at x = +10 (t = 110).
    enter, length = _beam.cone_path((-100.0, 0.0, 10.0), (1.0, 0.0, 0.0), Z_AXIS, 45.0)
    assert enter == pytest.approx(90.0, rel=1e-12)
    assert length == pytest.approx(20.0, rel=1e-12)


def test_beam_that_misses_cone_is_rejected():
    with pytest.raises(ValueError, match="does not cross"):
        _beam.cone_path((-100.0, 0.0, 10.0), (0.0, 1.0, 0.0), Z_AXIS, 45.0)


def test_beam_pointing_away_from_cone_is_rejected():
    with pytest.raises(ValueError, match="does not cross"):
        _beam.cone_path((0.0, 0.0, -10.0), DOWN, Z_AXIS, 10.0)


def test_beam_inside_cone_forever_is_unbounded():
    with pytest.raises(ValueError, match="unbounded"):
        _beam.cone_path((0.0, 0.0, 10.0), Z_AXIS, Z_AXIS, 10.0)


def test_beam_parallel_to_cone_generator_is_unbounded():
    # Direction (1, 0, 1)/sqrt(2) is a generator of the 45 deg cone: a = 0 (linear case).
    direction = (1 / math.sqrt(2), 0.0, 1 / math.sqrt(2))
    with pytest.raises(ValueError, match="unbounded"):
        _beam.cone_path((-10.0, 0.0, 0.0), direction, Z_AXIS, 45.0)


def test_t_max_clips_the_path():
    enter, length = _beam.cone_path((0.0, 0.0, 100.0), DOWN, Z_AXIS, 10.0, t_max=40.0)
    assert enter == pytest.approx(0.0)
    assert length == pytest.approx(40.0)


def test_unit_vector_convention():
    assert _beam.unit_vector(90.0, 0.0) == pytest.approx((1.0, 0.0, 0.0), abs=1e-12)
    assert _beam.unit_vector(90.0, 90.0) == pytest.approx((0.0, -1.0, 0.0), abs=1e-12)
    assert _beam.unit_vector(0.0, 45.0) == pytest.approx(Z_AXIS)
    with pytest.raises(ValueError, match="zenith"):
        _beam.unit_vector(200.0, 0.0)


# laser-beam


def test_laser_beam_is_one_line_with_cone_zero_and_beta_one():
    description = laser_beam.build(parse(laser_beam, LASER))
    assert len(description["events"]) == 1
    assert description["events"][0]["repeat"] == 1
    line = emission(description)
    assert line["kind"] == "line"
    assert line["beta"] == 1.0
    assert line["angular"] == {"kind": "cone", "angle_deg": 0.0}
    assert line["target"] == "none"
    assert line["spectrum"] == {"kind": "fixed", "wavelength_nm": 355.0}
    assert description["telescopes"] == [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 12.5}]


def test_laser_beam_enters_and_leaves_the_fov_cone():
    line = emission(laser_beam.build(parse(laser_beam, LASER)))
    tan_half = math.tan(math.radians(7.5))
    start = line["start_m"]
    end = [start[i] + line["direction"][i] * line["length_m"] for i in range(3)]
    assert on_cone_surface(start, tan_half)
    assert on_cone_surface(end, tan_half)
    assert line["length_m"] > 0.0


def test_laser_beam_photons_are_per_metre_times_length():
    line = emission(laser_beam.build(parse(laser_beam, LASER + " --photons-per-m 2.5e3")))
    assert line["photons"] == pytest.approx(2.5e3 * line["length_m"], rel=1e-12)


def test_laser_beam_along_axis_length_is_path_to_apex():
    text = (
        "--observation-level 2150 --spectrum fixed:wavelength_nm=355 "
        "--laser-position 0,0,100 --laser-zenith 180 --laser-azimuth 0 --photons-per-m 1e4 "
        "--telescope-zenith 0 --telescope-azimuth 0 --telescope-radius 12.5 --fov-deg 10 "
        "--events 3"
    )
    description = laser_beam.build(parse(laser_beam, text))
    assert description["events"][0]["repeat"] == 3
    line = emission(description)
    assert line["start_m"] == pytest.approx([0.0, 0.0, 100.0], abs=1e-9)
    assert line["length_m"] == pytest.approx(100.0)
    assert line["photons"] == pytest.approx(1e6)


def test_laser_beam_pointing_away_from_the_fov_is_rejected():
    text = LASER.replace("--laser-position 1000,0,0", "--laser-position 0,0,-100").replace(
        "--laser-zenith 30 --laser-azimuth 180", "--laser-zenith 180 --laser-azimuth 0"
    )
    with pytest.raises(ValueError, match="does not cross"):
        laser_beam.build(parse(laser_beam, text))


@pytest.mark.parametrize(
    "option",
    [
        "--observation-level",
        "--spectrum",
        "--laser-position",
        "--laser-zenith",
        "--laser-azimuth",
        "--photons-per-m",
        "--telescope-zenith",
        "--telescope-azimuth",
        "--telescope-radius",
        "--fov-deg",
    ],
)
def test_laser_beam_required_options_are_rejected(option):
    with pytest.raises(SystemExit):
        parse(laser_beam, without(LASER, option))


def test_laser_beam_fov_must_be_below_90():
    with pytest.raises(ValueError, match="field of view"):
        laser_beam.build(parse(laser_beam, LASER.replace("--fov-deg 15", "--fov-deg 90")))


# lidar


def test_photons_per_pulse_uses_exact_si_constants():
    expected = 1e-6 * 532e-9 / (6.62607015e-34 * 299792458.0)
    assert lidar.photons_per_pulse(1e-6, 532.0) == pytest.approx(expected, rel=1e-14)
    assert lidar.photons_per_pulse(1e-6, 532.0) == pytest.approx(2.6782e12, rel=1e-4)


def test_lidar_clipped_photons_scale_with_path_fraction():
    line = emission(lidar.build(parse(lidar, LIDAR)))
    # Beam straight down at x = 1000 m; cone of half angle 40 deg about z enters at
    # z = 1000 / tan(40 deg) and the satellite is at z = 4000 m. Full length: 4000 m.
    clipped = 4000.0 - 1000.0 / math.tan(math.radians(40.0))
    assert line["length_m"] == pytest.approx(clipped, rel=1e-12)
    assert line["start_m"] == pytest.approx([1000.0, 0.0, 4000.0], abs=1e-9)
    per_pulse = lidar.photons_per_pulse(1e-6, 532.0)
    assert line["photons"] == pytest.approx(per_pulse * clipped / 4000.0, rel=1e-12)


def test_lidar_structure_and_divergence():
    description = lidar.build(parse(lidar, LIDAR))
    line = emission(description)
    assert line["kind"] == "line"
    assert line["beta"] == 1.0
    assert line["target"] == "none"
    assert line["spectrum"] == {"kind": "fixed", "wavelength_nm": 532.0}
    assert line["angular"] == {"kind": "filled_cone", "angle_deg": pytest.approx(0.05)}
    assert line["direction"] == pytest.approx(DOWN, abs=1e-12)
    assert description["events"][0]["repeat"] == 1


def test_lidar_zero_divergence_is_cone_zero():
    line = emission(lidar.build(parse(lidar, LIDAR.replace("--divergence 0.1", "--divergence 0"))))
    assert line["angular"] == {"kind": "cone", "angle_deg": 0.0}


def test_lidar_on_axis_pulse_is_fully_inside_the_cone():
    text = (
        "--observation-level 2150 --satellite-altitude 400000 --beam-zenith 180 "
        "--beam-azimuth 0 --divergence 0 --pulse-energy 2e-3 --wavelength-nm 355 "
        "--telescope-zenith 0 --telescope-azimuth 0 --telescope-radius 12.5 --fov-deg 2 "
        "--events 5"
    )
    description = lidar.build(parse(lidar, text))
    assert description["events"][0]["repeat"] == 5
    line = emission(description)
    assert line["start_m"] == pytest.approx([0.0, 0.0, 400000.0], abs=1e-6)
    assert line["length_m"] == pytest.approx(400000.0)
    assert line["photons"] == pytest.approx(lidar.photons_per_pulse(2e-3, 355.0), rel=1e-12)


def test_lidar_beam_pointing_up_is_rejected():
    with pytest.raises(ValueError, match="towards the ground"):
        lidar.build(parse(lidar, LIDAR.replace("--beam-zenith 180", "--beam-zenith 60")))


def test_lidar_satellite_missing_the_fov_is_rejected():
    text = LIDAR.replace("--satellite-position 1000,0,4000", "--satellite-position 5000,0,4000")
    with pytest.raises(ValueError, match="does not cross"):
        lidar.build(parse(lidar, text.replace("--fov-deg 80", "--fov-deg 10")))


@pytest.mark.parametrize(
    "option",
    [
        "--observation-level",
        "--beam-zenith",
        "--beam-azimuth",
        "--divergence",
        "--pulse-energy",
        "--wavelength-nm",
        "--telescope-zenith",
        "--telescope-azimuth",
        "--telescope-radius",
        "--fov-deg",
    ],
)
def test_lidar_required_options_are_rejected(option):
    with pytest.raises(SystemExit):
        parse(lidar, without(LIDAR, option))


def test_lidar_satellite_position_or_altitude_exactly_one():
    with pytest.raises(SystemExit):
        parse(lidar, LIDAR + " --satellite-altitude 4000")
    with pytest.raises(SystemExit):
        parse(lidar, without(LIDAR, "--satellite-position"))


# night-sky


def test_night_sky_structure(tmp_path):
    table = tmp_path / "sky.dat"
    table.write_text("# wavelength density\n300 1\n600 1\n")
    description = night_sky.build(parse(night_sky, NIGHT_SKY.format(table=table)))
    assert description["telescopes"] == [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 12.5}]
    assert description["events"][0]["repeat"] == 1
    background = emission(description)
    assert background["kind"] == "sky_background"
    assert background["centre_m"] == [0.0, 0.0, 0.0]
    assert background["radius_m"] == 12.5
    assert background["axis"] == pytest.approx(DOWN)
    assert background["viewcone_deg"] == 10.0
    assert background["photons"] == 1e6
    assert background["spectrum"] == {
        "kind": "table",
        "wavelength_nm": [300.0, 600.0],
        "density": [1.0, 1.0],
    }


def test_night_sky_axis_is_minus_telescope_pointing(tmp_path):
    table = tmp_path / "sky.dat"
    table.write_text("300 1\n600 1\n")
    text = NIGHT_SKY.format(table=table) + " --telescope-zenith 30 --telescope-azimuth 0"
    background = emission(night_sky.build(parse(night_sky, text)))
    # Pointing (sin 30, 0, cos 30) towards north; the axis is its negative.
    assert background["axis"] == pytest.approx((-0.5, 0.0, -math.sqrt(3) / 2), abs=1e-12)


@pytest.mark.parametrize(
    "option",
    [
        "--observation-level",
        "--spectrum",
        "--photons",
        "--viewcone-deg",
        "--telescope-position",
        "--telescope-radius",
    ],
)
def test_night_sky_required_options_are_rejected(option, tmp_path):
    table = tmp_path / "sky.dat"
    table.write_text("300 1\n600 1\n")
    with pytest.raises(SystemExit):
        parse(night_sky, without(NIGHT_SKY.format(table=table), option))


def test_night_sky_rejects_invalid_view_cone(tmp_path):
    table = tmp_path / "sky.dat"
    table.write_text("300 1\n600 1\n")
    text = NIGHT_SKY.format(table=table).replace("--viewcone-deg 10", "--viewcone-deg 0")
    with pytest.raises(ValueError, match="view cone"):
        night_sky.build(parse(night_sky, text))


# Normalized run descriptions (only if the config module exists).


def test_normalize_accepts_preset_output(tmp_path):
    config = pytest.importorskip("obdeect.sources.config")
    table = tmp_path / "sky.dat"
    table.write_text("300 1\n600 1\n")
    for module, text in [
        (laser_beam, LASER),
        (lidar, LIDAR),
        (night_sky, NIGHT_SKY.format(table=table)),
    ]:
        normalized = config.normalize(module.build(parse(module, text)))
        assert normalized["events"][0]["emissions"]
