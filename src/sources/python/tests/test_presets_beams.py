"""Pulsed-laser presets (lidar, laser-beam) and night-sky: run descriptions from arguments."""

import argparse
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from obdeect.sources.presets import _pulsed_beam, laser_beam, lidar, night_sky

Z_AXIS = (0.0, 0.0, 1.0)
DOWN = (0.0, 0.0, -1.0)
C_M_NS = 0.299792458
H_C = 6.62607015e-34 * 299792458.0

LASER = (
    "--observation-level 2150 --laser-position 1000,0,4000 --beam-zenith 180 --beam-azimuth 0 "
    "--divergence 0.1 --pulse-energy 1e-6 --wavelength-nm 355 --telescope-radius 12.5"
)
LIDAR = (
    "--observation-level 2150 --satellite-altitude 400000 --satellite-offset 1000,0 "
    "--beam-zenith 180 --beam-azimuth 0 --divergence 0.1 --pulse-energy 1e-6 "
    "--wavelength-nm 532 --telescope-radius 12.5"
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


def replaced(text, option, value):
    """Option text with the value of one ``--option VALUE`` pair replaced."""
    tokens = text.split()
    tokens[tokens.index(option) + 1] = value
    return " ".join(tokens)


def emission(description):
    """The single emission of the first event of a run description."""
    emissions = description["events"][0]["emissions"]
    assert len(emissions) == 1
    return emissions[0]


def build_lidar(text):
    return lidar.build(parse(lidar, text))


def build_laser(text):
    return laser_beam.build(parse(laser_beam, text))


# Shared builder: direction convention, photons, bunches


def test_unit_vector_convention():
    assert _pulsed_beam.unit_vector(90.0, 0.0) == pytest.approx((1.0, 0.0, 0.0), abs=1e-12)
    assert _pulsed_beam.unit_vector(90.0, 90.0) == pytest.approx((0.0, -1.0, 0.0), abs=1e-12)
    assert _pulsed_beam.unit_vector(0.0, 45.0) == pytest.approx(Z_AXIS)
    assert _pulsed_beam.unit_vector(180.0, 0.0) == pytest.approx(DOWN)
    with pytest.raises(ValueError, match="zenith"):
        _pulsed_beam.unit_vector(200.0, 0.0)


def test_photon_number_is_exact_si_count():
    expected = 1e-6 * 532e-9 / H_C
    assert _pulsed_beam.photons_per_pulse(1e-6, 532.0) == pytest.approx(expected, rel=1e-14)
    assert emission(build_lidar(LIDAR))["photons"] == pytest.approx(expected, rel=1e-12)
    assert expected == pytest.approx(2.6782e12, rel=1e-4)


def test_photon_number_is_independent_of_rows():
    counts = {emission(build_lidar(LIDAR + f" --rows {rows}"))["photons"] for rows in (1e2, 1e5)}
    assert len(counts) == 1


def test_bunch_size_follows_rows():
    photons = emission(build_lidar(LIDAR))["photons"]
    assert emission(build_lidar(LIDAR))["bunch_size"] == pytest.approx(photons / 1e4, rel=1e-12)
    line = emission(build_lidar(LIDAR + " --rows 1000"))
    assert line["bunch_size"] == pytest.approx(photons / 1e3, rel=1e-12)
    assert emission(build_lidar(LIDAR + " --rows 1e30"))["bunch_size"] == 1.0


def test_laser_beam_is_one_downward_point_emission():
    description = build_laser(LASER)
    assert description["telescopes"] == [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 12.5}]
    assert description["events"][0]["repeat"] == 1
    point = emission(description)
    assert point["kind"] == "point"
    assert point["position_m"] == pytest.approx([1000.0, 0.0, 4000.0])
    assert point["axis"] == pytest.approx(DOWN, abs=1e-12)
    assert point["angular"] == {"kind": "filled_cone", "angle_deg": pytest.approx(0.05)}
    assert point["spectrum"] == {"kind": "fixed", "wavelength_nm": 355.0}
    assert point["pulse"] == {"kind": "fixed", "time_ns": 0.0}
    assert point["target"] == "none"
    assert point["time_offset_ns"] == 0.0
    assert point["to_plane_z0"] is False


def test_lidar_structure_and_zero_divergence():
    point = emission(build_lidar(LIDAR))
    assert point["kind"] == "point"
    assert point["position_m"] == pytest.approx([1000.0, 0.0, 400000.0])
    assert point["spectrum"] == {"kind": "fixed", "wavelength_nm": 532.0}
    point = emission(build_lidar(LIDAR.replace("--divergence 0.1", "--divergence 0")))
    assert point["angular"] == {"kind": "cone", "angle_deg": 0.0}


def test_pulse_profile_is_passed_through():
    point = emission(build_lidar(LIDAR + " --pulse gauss:sigma_ns=2.5"))
    assert point["pulse"] == {"kind": "gauss", "sigma_ns": 2.5}


def test_run_metadata_records_the_pulse_and_beam():
    metadata = build_lidar(LIDAR)["run"]["metadata"]
    assert metadata["pulse_energy_j"] == 1e-6
    assert metadata["wavelength_nm"] == 532.0
    assert metadata["divergence_deg"] == 0.1
    assert metadata["beam_zenith_deg"] == 180.0
    assert metadata["beam_azimuth_deg"] == 0.0
    assert metadata["laser_x_m"] == 1000.0
    assert metadata["laser_y_m"] == 0.0
    assert metadata["laser_z_m"] == 400000.0
    assert metadata["start_height_m"] == 400000.0


# Start height


def test_start_height_moves_the_point_along_the_beam():
    text = LIDAR.replace("--beam-zenith 180", "--beam-zenith 170").replace(
        "--beam-azimuth 0", "--beam-azimuth 30"
    )
    description = build_lidar(text + " --start-height 100000")
    point = emission(description)
    laser = [1000.0, 0.0, 400000.0]
    beam = _pulsed_beam.unit_vector(170.0, 30.0)
    assert point["axis"] == pytest.approx(list(beam), abs=1e-12)
    offset = [point["position_m"][i] - laser[i] for i in range(3)]
    path = math.sqrt(sum(o * o for o in offset))
    assert point["position_m"][2] == pytest.approx(100000.0, abs=1e-6)
    assert offset == pytest.approx([path * b for b in beam], rel=1e-9, abs=1e-6)
    assert point["time_offset_ns"] == pytest.approx(path / C_M_NS, rel=1e-12)
    assert description["run"]["metadata"]["start_height_m"] == pytest.approx(100000.0)


def test_start_height_zero_means_the_laser_itself():
    point = emission(build_lidar(LIDAR))
    assert point["time_offset_ns"] == 0.0


@pytest.mark.parametrize("height", ["400000", "500000", "-1"])
def test_start_height_outside_the_beam_is_rejected(height):
    with pytest.raises(ValueError, match="start-height"):
        build_lidar(LIDAR + f" --start-height {height}")


# Direction and tilt


@pytest.mark.parametrize("zenith", ["60", "90", "0"])
def test_upward_or_horizontal_beam_is_rejected(zenith):
    with pytest.raises(ValueError, match="point downwards"):
        build_lidar(replaced(LIDAR, "--beam-zenith", zenith))


def test_too_tilted_beam_is_rejected_with_reason():
    with pytest.raises(ValueError, match="plane-parallel"):
        build_lidar(replaced(LIDAR, "--beam-zenith", "100"))


def test_tilt_limit_is_inclusive_and_configurable():
    build_lidar(replaced(LIDAR, "--beam-zenith", "120"))  # tilt 60 deg, the default limit
    with pytest.raises(ValueError, match="max-tilt"):
        build_laser(LASER.replace("--beam-zenith 180", "--beam-zenith 115") + " --max-tilt 30")
    build_laser(LASER.replace("--beam-zenith 180", "--beam-zenith 115") + " --max-tilt 65")


def test_max_tilt_must_be_below_90():
    with pytest.raises(ValueError, match="max-tilt"):
        build_lidar(LIDAR + " --max-tilt 90")


def test_laser_below_the_observation_level_is_rejected():
    with pytest.raises(ValueError, match="above the observation level"):
        build_laser(LASER.replace("--laser-position 1000,0,4000", "--laser-position 1000,0,0"))


def test_divergence_must_be_below_180():
    with pytest.raises(ValueError, match="divergence"):
        build_lidar(replaced(LIDAR, "--divergence", "180"))


# Direct beam and telescopes


def test_direct_beam_hitting_a_telescope_is_rejected():
    text = LIDAR.replace("--satellite-offset 1000,0", "--satellite-offset 0,0")
    with pytest.raises(ValueError, match="the direct laser beam hits telescope 1"):
        build_lidar(text)


def test_direct_beam_hit_is_kept_with_allow_direct_beam():
    text = LIDAR.replace("--satellite-offset 1000,0", "--satellite-offset 0,0")
    description = build_lidar(text + " --allow-direct-beam")
    assert emission(description)["kind"] == "point"


def test_direct_beam_missing_a_telescope_is_accepted():
    assert emission(build_lidar(LIDAR))["kind"] == "point"  # 1000 m off axis


def test_divergence_cone_radius_counts_in_the_direct_beam_check():
    # Same 1000 m offset, but a 10 deg beam at 400 km has a radius of ~35 km at the telescope.
    with pytest.raises(ValueError, match="hits telescope 1"):
        build_lidar(LIDAR.replace("--divergence 0.1", "--divergence 10"))


def test_telescope_above_the_start_is_not_hit_by_the_beam():
    text = without(
        LIDAR.replace("--satellite-offset 1000,0", "--satellite-offset 0,0"), "--telescope-radius"
    )
    description = build_lidar(text + " --telescope 0,0,500000,12.5")
    assert description["telescopes"] == [
        {"id": 1, "position_m": [0.0, 0.0, 500000.0], "radius_m": 12.5}
    ]


def test_several_telescopes_get_ids_and_the_hit_is_named():
    text = without(
        LIDAR.replace("--satellite-offset 1000,0", "--satellite-offset 0,0"), "--telescope-radius"
    )
    with pytest.raises(ValueError, match="hits telescope 2"):
        build_lidar(text + " --telescope 1000,0,0,12.5 --telescope 0,0,0,12.5")
    description = build_lidar(
        text + " --telescope 1000,0,0,12.5 --telescope -1000,0,0,12.5 --telescope 0,0,0,1"
        " --allow-direct-beam"
    )
    assert [t["id"] for t in description["telescopes"]] == [1, 2, 3]
    assert description["telescopes"][1] == {
        "id": 2,
        "position_m": [-1000.0, 0.0, 0.0],
        "radius_m": 12.5,
    }


def test_telescope_option_needs_four_numbers():
    with pytest.raises(SystemExit):
        parse(lidar, LIDAR.replace("--telescope-radius 12.5", "--telescope 1,2,3"))


# Options


def test_lidar_satellite_offset_only_with_altitude():
    text = without(LIDAR, "--satellite-offset").replace(
        "--satellite-altitude 400000", "--satellite-position 1000,0,400000"
    )
    assert emission(build_lidar(text))["position_m"] == pytest.approx([1000.0, 0.0, 400000.0])
    with pytest.raises(ValueError, match="satellite-offset"):
        build_lidar(text + " --satellite-offset 1,2")


def test_lidar_satellite_position_or_altitude_exactly_one():
    with pytest.raises(SystemExit):
        parse(lidar, LIDAR + " --satellite-position 1,2,3")
    with pytest.raises(SystemExit):
        parse(
            lidar, without(LIDAR, "--satellite-altitude").replace("--satellite-offset 1000,0", "")
        )


@pytest.mark.parametrize(
    "option",
    [
        "--observation-level",
        "--beam-zenith",
        "--beam-azimuth",
        "--divergence",
        "--pulse-energy",
        "--wavelength-nm",
        "--telescope-radius",
    ],
)
def test_lidar_required_options_are_rejected(option):
    with pytest.raises(SystemExit):
        parse(lidar, without(LIDAR, option))


@pytest.mark.parametrize(
    "option",
    [
        "--observation-level",
        "--laser-position",
        "--beam-zenith",
        "--beam-azimuth",
        "--divergence",
        "--pulse-energy",
        "--wavelength-nm",
        "--telescope-radius",
    ],
)
def test_laser_beam_required_options_are_rejected(option):
    with pytest.raises(SystemExit):
        parse(laser_beam, without(LASER, option))


def test_rows_must_be_positive():
    with pytest.raises(ValueError, match="rows"):
        build_lidar(LIDAR + " --rows 0")


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


# Normalized run descriptions and the end-to-end chain


@pytest.mark.parametrize("module, text", [(laser_beam, LASER), (lidar, LIDAR)])
def test_normalize_accepts_preset_output(module, text):
    config = pytest.importorskip("obdeect.sources.config")
    normalized = config.normalize(module.build(parse(module, text)))
    assert normalized["events"][0]["emissions"][0]["kind"] == "point"


def _obdeect_command():
    """The obdeect command of this interpreter's environment, or None."""
    candidate = Path(sys.executable).with_name("obdeect")
    if candidate.exists():
        return str(candidate)
    return shutil.which("obdeect")


def _run(*args):
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    return result.stdout + result.stderr


@pytest.fixture
def obdeect_cli():
    command = _obdeect_command()
    if command is None:
        pytest.skip("obdeect command not available")
    probe = subprocess.run(
        [command, "propagation", "run", "--help"], capture_output=True, check=False
    )
    if probe.returncode != 0:
        pytest.skip("obdeect propagation command not available")
    return command


def _photon_rows(command, path):
    """Photons of a stream file, from ``obdeect stream cat --summary``."""
    match = re.search(
        r"^total: (\d+) photons", _run(command, "stream", "cat", "--summary", path), re.MULTILINE
    )
    assert match is not None
    return int(match.group(1))


def test_lidar_chain_to_scattered_light(obdeect_cli, tmp_path):
    data = Path(__file__).resolve().parents[3] / "propagation" / "cpp" / "tests" / "data"
    primaries = tmp_path / "primaries.pst"
    scattered = tmp_path / "scattered.pst"
    _run(
        obdeect_cli, "sources", "run", "lidar",
        "--observation-level", "2000", "--satellite-altitude", "100000",
        "--satellite-offset", "1000,0", "--beam-zenith", "180", "--beam-azimuth", "0",
        "--divergence", "0.01", "--pulse-energy", "1e-3", "--wavelength-nm", "355",
        "--telescope-radius", "12.5", "--rows", "1000", "-o", str(primaries),
    )  # fmt: skip
    _run(
        obdeect_cli, "propagation", "run",
        "--transmission", str(data / "atm_trans_2000_1_10_0_0_2000.dat"),
        "--site-altitude", "2000", "--seed", "1", "--mode", "weighted",
        "--scattering", str(data / "scattering_exponential.dat"), "--scattering-points", "64",
        "-i", str(primaries), "-o", str(scattered),
    )  # fmt: skip
    assert _photon_rows(obdeect_cli, scattered) > _photon_rows(obdeect_cli, primaries)
