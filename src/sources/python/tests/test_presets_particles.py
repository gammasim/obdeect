"""Cherenkov helpers and the muon and shower presets: formulas, structure and option errors."""

import argparse
import math

import numpy as np
import pytest
from obdeect.sources.presets import _cherenkov, muon, shower

C_M_PER_NS = 0.299792458
ALPHA = 0.0072973526

MUON_ARGV = [
    "--observation-level", "2150",
    "--telescope-radius", "12.5",
    "--start-height", "3000",
    "--n-minus-one", "0.000292",
    "--wavelength-min", "300",
    "--wavelength-max", "600",
    "--events", "3",
    "--seed", "7",
]  # fmt: skip

SHOWER_ARGV = [
    "--observation-level", "2150",
    "--telescope-radius", "12.5",
    "--wavelength-min", "300",
    "--wavelength-max", "600",
    "--n-minus-one", "0.0003",
    "--start-height", "30000",
    "--zenith", "0",
    "--azimuth", "0",
    "--step", "500",
    "--nmax", "1e6",
    "--xmax", "600",
    "--x0", "100",
    "--gh-lambda", "70",
    "--atmosphere-depth0", "1030",
    "--scale-height", "8400",
    "--fluorescence-yield", "4",
    "--fluorescence-spectrum", "fixed:wavelength_nm=337",
    "--events", "2",
]  # fmt: skip


def parse(module, argv):
    parser = argparse.ArgumentParser()
    module.add_arguments(parser)
    return parser.parse_args(argv)


def frank_tamm_analytic(cos_theta, l1, l2, charge=1.0):
    return 2.0 * math.pi * ALPHA * charge**2 * (1.0 - cos_theta**2) * (1 / l1 - 1 / l2) * 1e9


def test_frank_tamm_hand_value():
    # 2 pi alpha = 0.0458506; x sin^2 = 0.64; x (1/300 - 1/600) nm^-1 = 1/600 nm^-1 = 1.6667e6 /m
    assert _cherenkov.frank_tamm_yield(0.6, 300.0, 600.0) == pytest.approx(48907.3, rel=1e-5)
    assert _cherenkov.frank_tamm_yield(0.6, 300.0, 600.0, charge=2.0) == pytest.approx(
        4 * 48907.3, rel=1e-5
    )


def test_frank_tamm_rejects_bad_input():
    with pytest.raises(ValueError):
        _cherenkov.frank_tamm_yield(1.01, 300.0, 600.0)
    with pytest.raises(ValueError):
        _cherenkov.frank_tamm_yield(0.6, 600.0, 300.0)


def test_index_profiles():
    assert _cherenkov.cherenkov_cos(0.000292, 1.0) == pytest.approx(1 / 1.000292)
    assert _cherenkov.cherenkov_cos(0.000292, 0.5) == pytest.approx(2 / 1.000292)
    exponential = _cherenkov.exponential_index(3e-4, 8000.0)
    assert exponential(0.0) == pytest.approx(3e-4)
    assert exponential(8000.0) == pytest.approx(3e-4 / math.e)
    table = _cherenkov.table_index([0.0, 1000.0, 3000.0], [3e-4, 2e-4, 1e-4])
    assert table(500.0) == pytest.approx(2.5e-4)
    assert table(-10.0) == pytest.approx(3e-4)
    assert table(1e5) == pytest.approx(1e-4)
    with pytest.raises(ValueError):
        _cherenkov.table_index([0.0, 0.0], [3e-4, 2e-4])
    with pytest.raises(ValueError):
        _cherenkov.constant_index(0.0)


def test_track_segments_split():
    segments = _cherenkov.track_segments([0.0, 0.0, 0.0], [0.0, 0.0, -1.0], 100.0, 30.0)
    assert len(segments) == 4
    assert all(s.length_m == pytest.approx(25.0) for s in segments)
    assert [s.distance_m for s in segments] == pytest.approx([25.0, 50.0, 75.0, 100.0])
    assert segments[-1].start_m == pytest.approx((0.0, 0.0, 100.0))
    assert segments[0].altitude_m == pytest.approx(12.5)


def test_time_convention_arrival_at_end_is_zero():
    assert _cherenkov.time_offset_ns(0.0, 1.0) == 0.0
    assert _cherenkov.time_offset_ns(3000.0, 1.0) == pytest.approx(-3000.0 / C_M_PER_NS)
    assert _cherenkov.time_offset_ns(3000.0, 0.5) == pytest.approx(-6000.0 / C_M_PER_NS)


def test_muon_build_structure():
    description = muon.build(parse(muon, MUON_ARGV))
    assert set(description) == {"run", "telescopes", "events"}
    assert description["run"]["observation_level_m"] == 2150.0
    assert description["run"]["frame"] == "sky"
    assert description["telescopes"] == [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": 12.5}]
    assert len(description["events"]) == 1
    event = description["events"][0]
    assert event["repeat"] == 3
    assert event["metadata"]["primary_id"] == 6
    assert len(event["emissions"]) == 1
    line = event["emissions"][0]
    assert line["kind"] == "line"
    assert line["start_m"] == pytest.approx([0.0, 0.0, 3000.0])
    assert line["direction"] == pytest.approx([0.0, 0.0, -1.0])
    assert line["length_m"] == pytest.approx(3000.0)
    assert line["beta"] == 1.0
    assert line["time_offset_ns"] == pytest.approx(-3000.0 / C_M_PER_NS)
    assert line["angular"]["kind"] == "cone"
    assert line["angular"]["angle_deg"] == pytest.approx(math.degrees(math.acos(1 / 1.000292)))
    assert line["spectrum"] == {"kind": "inverse_square", "low_nm": 300.0, "high_nm": 600.0}
    assert line["target"] == "auto"


def test_muon_total_photons_analytic():
    description = muon.build(parse(muon, MUON_ARGV))
    cos_theta = 1 / 1.000292
    expected_per_event = frank_tamm_analytic(cos_theta, 300.0, 600.0) * 3000.0
    total = sum(e["photons"] for e in description["events"][0]["emissions"])
    assert total == pytest.approx(expected_per_event, rel=1e-9)
    assert total * 3 == pytest.approx(expected_per_event * 3, rel=1e-9)


def test_muon_zenith_length_and_profile_integral(tmp_path):
    # Zenith 30 deg: track length h / cos(z). Index table: photons = integral of the yield.
    table = tmp_path / "index.txt"
    table.write_text("0 3e-4\n1000 2e-4\n2000 1e-4\n3000 5e-5\n10000 1e-5\n")
    argv = [a for a in MUON_ARGV + ["--zenith", "30"] if a not in ("--n-minus-one", "0.000292")]
    argv += ["--index-profile", str(table)]
    description = muon.build(parse(muon, argv))
    emissions = description["events"][0]["emissions"]
    # Path length 3000 / cos(30 deg) = 3464 m, in ceil(3464 / 10) = 347 steps.
    assert len(emissions) == 347
    assert sum(e["length_m"] for e in emissions) == pytest.approx(
        3000.0 / math.cos(math.radians(30))
    )
    total = sum(e["photons"] for e in emissions)

    # Fine midpoint integral along the track, with the same piecewise linear index table.
    h = np.linspace(0.0, 3000.0, 3_000_001)
    h_mid = 0.5 * (h[1:] + h[:-1])
    eta = np.interp(h_mid, [0, 1000, 2000, 3000, 10000], [3e-4, 2e-4, 1e-4, 5e-5, 1e-5])
    cos_theta = 1.0 / (1.0 + eta)
    density = 2.0 * math.pi * ALPHA * (1.0 - cos_theta**2) * (1 / 300.0 - 1 / 600.0) * 1e9
    integral = float(np.sum(density) * (h[1] - h[0])) / math.cos(math.radians(30.0))
    assert total == pytest.approx(integral, rel=1e-3)


def test_muon_time_convention_profile(tmp_path):
    table = tmp_path / "index.txt"
    table.write_text("0 3e-4\n10000 1e-5\n")
    argv = [a for a in MUON_ARGV if a not in ("--n-minus-one", "0.000292")]
    argv += ["--index-profile", str(table)]
    emissions = muon.build(parse(muon, argv))["events"][0]["emissions"]
    offsets = [e["time_offset_ns"] for e in emissions]
    assert offsets[0] == pytest.approx(-10.0 / C_M_PER_NS)  # lowest step: 10 m from the impact
    assert offsets[-1] == pytest.approx(-3000.0 / C_M_PER_NS)  # top step: full track length
    assert offsets == sorted(offsets, reverse=True)


def test_muon_core_scatter_events_are_separate_and_deterministic():
    argv = MUON_ARGV + ["--core-scatter-radius", "100", "--core-scatter-count", "5"]
    first = muon.build(parse(muon, argv))
    second = muon.build(parse(muon, argv))
    assert len(first["events"]) == 5
    cores = [(e["metadata"]["core_x_m"], e["metadata"]["core_y_m"]) for e in first["events"]]
    assert all(e["repeat"] == 3 for e in first["events"])
    assert all(math.hypot(x, y) <= 100.0 for x, y in cores)
    assert cores == [
        (e["metadata"]["core_x_m"], e["metadata"]["core_y_m"]) for e in second["events"]
    ]
    for event, (x, y) in zip(first["events"], cores):
        assert event["emissions"][0]["start_m"][:2] == pytest.approx([x, y])


def test_muon_option_errors():
    with pytest.raises(SystemExit):
        parse(muon, [a for a in MUON_ARGV if a not in ("--telescope-radius", "12.5")])
    with pytest.raises(SystemExit):
        parse(muon, [a for a in MUON_ARGV if a not in ("--start-height", "3000")])
    with pytest.raises(SystemExit):  # index: exactly one option
        parse(muon, MUON_ARGV + ["--index-exponential", "3e-4", "8000"])
    with pytest.raises(ValueError, match="zenith"):
        muon.build(parse(muon, MUON_ARGV + ["--zenith", "95"]))
    with pytest.raises(ValueError, match="unknown primary name"):
        muon.build(parse(muon, MUON_ARGV + ["--primary", "neutron"]))
    with pytest.raises(ValueError):
        muon.build(parse(muon, MUON_ARGV + ["--n-minus-one", "0"]))


def test_gaisser_hillas_peaks_at_xmax():
    xmax = 600.0
    values = [shower.gaisser_hillas(x, 1e6, xmax, 100.0, 70.0) for x in np.arange(101, 2000, 1.0)]
    assert shower.gaisser_hillas(xmax, 1e6, xmax, 100.0, 70.0) == pytest.approx(1e6)
    assert shower.gaisser_hillas(100.0, 1e6, xmax, 100.0, 70.0) == 0.0
    peak_x = 101.0 + int(np.argmax(values))
    assert abs(peak_x - xmax) <= 1.0


def test_shower_profile_peaks_at_xmax():
    description = shower.build(parse(shower, SHOWER_ARGV))
    emissions = description["events"][0]["emissions"]
    cherenkov = [e for e in emissions if e["angular"]["kind"] == "cone"]
    # Midpoint altitudes of the steps: top step at 30000 m, 500 m long, zenith 0.
    mids = [e["start_m"][2] - 0.5 * e["length_m"] for e in cherenkov]
    photons = [e["photons"] for e in cherenkov]
    h_xmax = 8400.0 * math.log(1030.0 / 600.0)  # X_v(h) = Xmax
    h_x0 = 8400.0 * math.log(1030.0 / 100.0)  # X_v(h) = X0: no light above
    assert max(mids) <= h_x0 + 1e-6
    peak_mid = mids[int(np.argmax(photons))]
    assert abs(peak_mid - h_xmax) <= 500.0


def test_fluorescence_lines_are_isotropic():
    description = shower.build(parse(shower, SHOWER_ARGV))
    emissions = description["events"][0]["emissions"]
    fluorescence = [e for e in emissions if e["angular"]["kind"] == "isotropic"]
    assert fluorescence
    for line in fluorescence:
        assert line["kind"] == "line"
        assert line["spectrum"] == {"kind": "fixed", "wavelength_nm": 337.0}
        assert line["beta"] == 1.0
        assert line["target"] == "auto"
        mid = line["start_m"][2] - 0.5 * line["length_m"]
        depth = shower.slant_depth(mid, 1030.0, 8400.0, 0.0)
        particles = shower.gaisser_hillas(depth, 1e6, 600.0, 100.0, 70.0)
        assert line["photons"] == pytest.approx(4.0 * line["length_m"] * particles, rel=1e-9)
    # Steps above X0 carry no particles, so the list starts at the lowest step (500 m from core).
    assert fluorescence[0]["time_offset_ns"] == pytest.approx(-500.0 / C_M_PER_NS)


def test_shower_option_errors():
    for missing in ("--fluorescence-yield", "--gh-lambda", "--step", "--zenith"):
        argv = list(SHOWER_ARGV)
        index = argv.index(missing)
        del argv[index : index + 2]
        with pytest.raises(SystemExit):
            parse(shower, argv)
    x0_argv = list(SHOWER_ARGV)
    x0_argv[x0_argv.index("--x0") + 1] = "700"
    with pytest.raises(ValueError, match="X0 < Xmax"):
        shower.build(parse(shower, x0_argv))
    bad = [a if a != "fixed:wavelength_nm=337" else "inverse_square:low_nm=300,high_nm=600"
           for a in SHOWER_ARGV]  # fmt: skip
    with pytest.raises(ValueError, match="fixed or table"):
        shower.build(parse(shower, bad))


def test_normalize_round_trip_if_available():
    config = pytest.importorskip("obdeect.sources.config")
    normalize = getattr(config, "normalize", None)
    if normalize is None:
        pytest.skip("obdeect.sources.config has no normalize")
    for module, argv in ((muon, MUON_ARGV), (shower, SHOWER_ARGV)):
        normalized = normalize(module.build(parse(module, argv)))
        assert normalized["events"]
