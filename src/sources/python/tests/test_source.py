"""Single emissions: kinds, angular, spectral and pulse profiles, chunks and validation."""

import math

import numpy as np
import pytest
from obdeect.sources.source import (
    CHUNK_ROWS,
    COLUMNS,
    emit,
    emit_chunk,
    expected_photons,
    expected_weight,
    iter_chunks,
    normalize_emission,
    photon_count,
)

FLASHER = {
    "kind": "point",
    "photons": 4e5,
    "position_m": [0.0, 0.0, 1.5],
    "axis": [0.0, 0.0, -1.0],
    "angular": {"kind": "filled_cone", "angle_deg": 11.0},
    "spectrum": {"kind": "fixed", "wavelength_nm": 385.0},
    "pulse": {"kind": "tophat", "width_ns": 2.5},
    "target": {
        "kind": "disc",
        "centre_m": [0.0, 0.0, 0.0],
        "normal": [0.0, 0.0, 1.0],
        "radius_m": 1.2,
    },
}

PLANE_WAVE = {
    "kind": "plane_wave",
    "photons": 1e5,
    "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
    "towards_source": [0.01, 0.0, 1.0],
    "disc_radius_m": 5.0,
    "disc_z_m": 0.0,
    "start_z_m": 30.0,
}


def point_with(angular, photons=2e5, **extra):
    spec = {
        "kind": "point",
        "photons": photons,
        "position_m": [0.0, 0.0, 10.0],
        "axis": [0.0, 0.0, -1.0],
        "angular": angular,
        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
    }
    spec.update(extra)
    return spec


def cosine_to_axis(chunk, axis):
    axis = np.asarray(axis, dtype=float)
    return (chunk["cx"] * axis[0] + chunk["cy"] * axis[1] + chunk["cz"] * axis[2]).astype(float)


def test_normalize_emission_fills_defaults():
    spec = {
        "kind": "point",
        "photons": 10,
        "position_m": [0, 0, 1],
        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
    }
    emission = normalize_emission(spec)
    assert emission["bunch_size"] == 1.0
    assert emission["target"] == {"kind": "none"}
    assert emission["pulse"] == {"kind": "fixed", "time_ns": 0.0}
    assert emission["photons"] == 10.0


def test_columns_are_float32_and_directions_unit():
    chunk = emit_chunk(FLASHER, seed=3, rows=5000)
    assert tuple(chunk) == COLUMNS
    assert all(array.dtype == np.float32 for array in chunk.values())
    length = np.sqrt(chunk["cx"] ** 2 + chunk["cy"] ** 2 + chunk["cz"] ** 2)
    assert np.allclose(length, 1.0, atol=1e-6)
    assert np.all((chunk["w"] >= 0.0) & (chunk["w"] <= 1.0))
    assert np.all(chunk["z"] == np.float32(1.5))
    assert np.all(chunk["wavelength"] == np.float32(385.0))


def test_disc_target_weights_sum_to_photons_reaching_the_disc():
    # The 11 deg cone (pointing at the disc) lands entirely on the 1.2 m disc: the weights
    # sum to the emitted number of photons (within the rounding of the count).
    photons = emit(FLASHER, seed=3)
    assert photons["w"].sum() == pytest.approx(4e5, rel=1e-3)
    assert len(photons["w"]) == photon_count(FLASHER, seed=3)
    assert expected_weight(FLASHER) == pytest.approx(4e5, rel=1e-3)


def test_plane_wave_starts_on_the_start_plane():
    # Infinite distance: rays start at z = start_z and are parallel to -towards_source; a
    # quarter of the disc (radius 2.5 of 5 m) is hit by 25 % of them.
    photons = emit(dict(PLANE_WAVE, disc_radius_m=5.0), seed=2)
    assert np.allclose(photons["z"], 30.0)
    towards = np.array([0.01, 0.0, 1.0]) / np.hypot(0.01, 1.0)
    assert np.allclose(photons["cx"], -towards[0], atol=1e-6)
    assert np.allclose(photons["cz"], -towards[2], atol=1e-6)
    s = 30.0 / photons["cz"].astype(float)
    hit = np.hypot(photons["x"] - s * photons["cx"], photons["y"] - s * photons["cy"])
    assert np.all(hit <= 5.0 + 1e-3)
    inside = np.count_nonzero(hit < 2.5) / len(hit)
    assert inside == pytest.approx(0.25, abs=4.0 * math.sqrt(0.25 * 0.75 / len(hit)))


def test_plane_wave_at_finite_distance_emits_from_its_position():
    spec = {
        "kind": "plane_wave",
        "photons": 20.0,
        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
        "towards_source": [0.6, 0.0, 0.8],
        "distance_m": 20.0,
        "disc_radius_m": 1.0,
        "disc_z_m": 1.5,
        "start_z_m": 30.0,
    }
    photons = emit(spec, seed=1)
    assert np.allclose(photons["x"], 12.0, atol=1e-4)
    assert np.allclose(photons["z"], 17.5, atol=1e-4)
    # t = (|start - source| - distance) / c with start = source: -distance / c.
    assert np.allclose(photons["t"], -20.0 / 0.299792458, atol=1e-3)


def test_line_emission_positions_lie_on_the_segment_and_cone_angle_holds():
    # Cherenkov-like cone of 1.5 deg about a vertical line from 100 km down to 90 km.
    spec = {
        "kind": "line",
        "photons": 2e4,
        "start_m": [0.0, 0.0, 1e5],
        "direction": [0.0, 0.0, -1.0],
        "length_m": 1e4,
        "angular": {"kind": "cone", "angle_deg": 1.5},
        "spectrum": {"kind": "fixed", "wavelength_nm": 532.0},
    }
    photons = emit(spec, seed=6)
    assert np.all((photons["z"] <= 1e5) & (photons["z"] >= 9e4))
    assert np.all(photons["x"] == 0.0)
    cosine = cosine_to_axis(photons, [0.0, 0.0, -1.0])
    assert np.allclose(cosine, math.cos(math.radians(1.5)), atol=1e-5)
    assert len(photons["w"]) == photon_count(spec, seed=6)


def test_line_profile_concentrates_photons_where_the_density_is():
    spec = {
        "kind": "line",
        "photons": 1e5,
        "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
        "start_m": [0.0, 0.0, 0.0],
        "direction": [1.0, 0.0, 0.0],
        "length_m": 10.0,
        "angular": {"kind": "isotropic"},
        "profile": {"distance_m": [0.0, 5.0, 10.0], "density": [0.0, 1.0, 0.0]},
    }
    photons = emit(spec, seed=4)
    assert np.all((photons["x"] >= 0.0) & (photons["x"] <= 10.0))
    assert photons["x"].mean() == pytest.approx(5.0, abs=0.05)


@pytest.mark.parametrize(
    ("angular", "expected_cosine"),
    [
        ({"kind": "cone", "angle_deg": 10.0}, math.cos(math.radians(10.0))),
        ({"kind": "lambertian"}, 2.0 / 3.0),  # E[cos theta] for intensity cos theta.
        ({"kind": "isotropic"}, 0.0),
    ],
)
def test_angular_profiles_mean_cosine(angular, expected_cosine):
    photons = emit(point_with(angular), seed=5)
    cosine = cosine_to_axis(photons, [0.0, 0.0, -1.0])
    sigma = max(np.std(cosine), 1e-12) / math.sqrt(len(cosine))
    assert cosine.mean() == pytest.approx(expected_cosine, abs=4.0 * sigma + 1e-6)


def test_gauss_pulse_and_exponential_moments():
    gauss = emit(
        point_with({"kind": "isotropic"}, pulse={"kind": "gauss", "sigma_ns": 2.0, "mean_ns": 5.0}),
        seed=4,
    )
    assert gauss["t"].mean() == pytest.approx(5.0, abs=4.0 * 2.0 / math.sqrt(len(gauss["t"])))
    assert gauss["t"].std() == pytest.approx(2.0, rel=0.02)
    decay = emit(
        point_with(
            {"kind": "isotropic"}, pulse={"kind": "exponential", "decay_ns": 3.0, "mean_ns": 1.0}
        ),
        seed=4,
    )
    assert decay["t"].mean() == pytest.approx(1.0, abs=4.0 * 3.0 / math.sqrt(len(decay["t"])))


def test_tabulated_spectrum_and_pulse_follow_their_tables():
    spec = point_with(
        {"kind": "isotropic"},
        spectrum={
            "kind": "table",
            "wavelength_nm": [300.0, 400.0, 500.0],
            "density": [0.0, 1.0, 0.0],
        },
        pulse={"kind": "table", "time_ns": [0.0, 4.0], "rate": [1.0, 1.0]},
    )
    photons = emit(spec, seed=9)
    assert np.all((photons["wavelength"] > 300.0) & (photons["wavelength"] < 500.0))
    assert photons["wavelength"].mean() == pytest.approx(400.0, abs=2.0)
    assert np.all((photons["t"] >= 0.0) & (photons["t"] <= 4.0))


def test_gauss_angular_table_matches_its_width():
    photons = emit(point_with({"kind": "gauss", "sigma_deg": 5.0}, photons=1e5), seed=2)
    angle = np.degrees(np.arccos(np.clip(cosine_to_axis(photons, [0.0, 0.0, -1.0]), -1, 1)))
    # Two-dimensional r.m.s. of the Gaussian angular profile: sigma.
    assert np.sqrt(np.mean(angle**2) / 2.0) == pytest.approx(5.0, rel=0.02)


def test_bunch_size_divides_the_candidates_and_multiplies_the_weights():
    plain = dict(FLASHER, photons=4e5)
    bunched = dict(FLASHER, photons=4e5, bunch_size=4.0)
    assert expected_photons(bunched) == pytest.approx(expected_photons(plain) / 4.0, rel=1e-12)
    assert expected_weight(bunched) == pytest.approx(expected_weight(plain), rel=1e-12)
    # A quarter of the candidates: the sampled weight sum scatters by a few per mille.
    photons = emit(bunched, seed=3)
    assert photons["w"].sum() == pytest.approx(4e5, rel=1e-2)


def test_time_offset_shifts_the_emission_time():
    base = emit(point_with({"kind": "isotropic"}, photons=1000.0), seed=1)
    shifted = emit(point_with({"kind": "isotropic"}, photons=1000.0, time_offset_ns=7.5), seed=1)
    assert np.allclose(shifted["t"] - base["t"], np.float32(7.5), atol=1e-4)


def test_chunks_do_not_depend_on_how_the_photons_are_split():
    whole = emit_chunk(FLASHER, seed=7, rows=CHUNK_ROWS)
    pieces = [
        emit_chunk(FLASHER, seed=7, first_index=first, rows=rows)
        for first, rows in ((0, 1000), (1000, 1), (1001, CHUNK_ROWS - 1001))
    ]
    for name in COLUMNS:
        joined = np.concatenate([piece[name] for piece in pieces])
        assert np.array_equal(joined, whole[name]), name
    again = emit(FLASHER, seed=7)
    assert np.array_equal(again["x"], emit(FLASHER, seed=7)["x"])


def test_iter_chunks_covers_all_photons_in_bounded_pieces():
    count = photon_count(FLASHER, seed=11)
    sizes = [len(chunk["w"]) for chunk in iter_chunks(FLASHER, seed=11)]
    assert sum(sizes) == count
    assert max(sizes) <= CHUNK_ROWS
    assert len(sizes) == math.ceil(count / CHUNK_ROWS)


def test_photon_count_is_the_expected_count_rounded():
    spec = dict(PLANE_WAVE, photons=2.5)
    counts = [photon_count(spec, seed=seed) for seed in range(4000)]
    assert set(counts) <= {2, 3}
    assert np.mean(counts) == pytest.approx(2.5, abs=4.0 * 0.5 / math.sqrt(len(counts)))
    # Candidates per photon of the 11 deg cone on the disc:
    # I_max / (2 pi integral(I du)) * pi R^2 / h^2 with integral(I du) = 1 - cos(11 deg).
    per_photon = 1.0 / (2.0 * math.pi * (1.0 - math.cos(math.radians(11.0))))
    per_photon *= math.pi * 1.2**2 / 1.5**2
    assert expected_photons(FLASHER) == pytest.approx(4e5 * per_photon, rel=1e-12)


def test_sphere_target_photons_all_hit_the_sphere():
    spec = point_with(
        {"kind": "isotropic"},
        photons=1e6,
        target={"kind": "sphere", "centre_m": [0.0, 0.0, 0.0], "radius_m": 2.0},
    )
    photons = emit(spec, seed=8)
    origin = np.stack([photons["x"], photons["y"], photons["z"]], axis=1).astype(float)
    direction = np.stack([photons["cx"], photons["cy"], photons["cz"]], axis=1).astype(float)
    along = -np.sum(origin * direction, axis=1)
    perpendicular = np.linalg.norm(-origin - direction * along[:, None], axis=1)
    assert np.all(perpendicular <= 2.0 + 1e-4)
    cos_max = math.sqrt(1.0 - 4.0 / 100.0)
    assert expected_photons(spec) == pytest.approx(1e6 * 0.5 * (1.0 - cos_max))


def test_auto_target_is_rejected_for_a_single_emission():
    with pytest.raises(ValueError, match="auto"):
        emit_chunk(dict(FLASHER, target={"kind": "auto"}), seed=1, rows=10)


def test_invalid_emissions_are_rejected():
    def rejected(spec, match=None):
        with pytest.raises(ValueError, match=match):
            emit_chunk(spec, seed=1, rows=10)

    rejected(dict(FLASHER, unknown=1), "unknown key")
    rejected({k: v for k, v in FLASHER.items() if k != "photons"}, "photons")
    rejected(dict(FLASHER, photons=-1.0))
    rejected(dict(FLASHER, spectrum={"kind": "nonsense"}), "unknown kind")
    rejected(dict(FLASHER, pulse={"kind": "gauss", "sigma_ns": 1.0, "delay": 2.0}), "unknown key")
    rejected(dict(FLASHER, pulse={"kind": "gauss", "mean_ns": 1.0}), "sigma_ns")
    rejected(
        point_with(
            {"kind": "cone", "angle_deg": 5.0},
            target={"kind": "disc", "centre_m": [0, 0, 0], "normal": [0, 0, 1], "radius_m": 1.0},
        ),
        "fixed-angle",
    )
    rejected(
        dict(
            PLANE_WAVE,
            target={"kind": "disc", "centre_m": [0, 0, 0], "normal": [0, 0, 1], "radius_m": 1.0},
        ),
        "no target",
    )
    rejected(dict(PLANE_WAVE, towards_source=[0.0, 0.0, 0.0]), "non-zero")
    with pytest.raises(ValueError, match="rows"):
        emit_chunk(FLASHER, seed=1, rows=CHUNK_ROWS + 1)
    with pytest.raises(ValueError, match="first_index"):
        emit_chunk(FLASHER, seed=1, first_index=-1)


@pytest.mark.parametrize("angle", [-10.0, 190.0, math.nan, math.inf])
def test_cone_angles_are_degrees_in_0_180(angle):
    # Negative and > 180 angles would alias to valid cones through 1 - cos.
    for kind in ("cone", "filled_cone"):
        with pytest.raises(ValueError, match="finite|in \\[0, 180\\]"):
            emit_chunk(point_with({"kind": kind, "angle_deg": angle}), seed=1, rows=10)
    assert expected_photons(point_with({"kind": "cone", "angle_deg": 180.0})) > 0.0


@pytest.mark.parametrize("angles", [[0.0, -5.0, 10.0], [0.0, 200.0], [0.0, math.nan, 30.0]])
def test_angular_table_angles_are_checked(angles):
    table = {"kind": "table", "angle_deg": angles, "value": [1.0] * len(angles)}
    with pytest.raises(ValueError):
        emit_chunk(point_with(table), seed=1, rows=10)
    ok = {"kind": "table", "angle_deg": [0.0, 90.0, 180.0], "value": [1.0, 0.5, 0.1]}
    assert expected_photons(point_with(ok)) > 0.0


def test_gauss_angular_width_must_be_representable():
    with pytest.raises(ValueError):
        emit_chunk(point_with({"kind": "gauss", "sigma_deg": 1e-200}), seed=1, rows=10)


def test_photon_counts_must_be_finite_and_exact():
    # Counts beyond 2^53 are not exact in a double; iter_chunks would loop without end.
    with pytest.raises(ValueError):
        photon_count(dict(FLASHER, photons=1e20), seed=1)
    with pytest.raises(ValueError):
        photon_count(dict(FLASHER, photons=math.inf), seed=1)
    # A disc target at h = 1e-170 makes the candidate count overflow.
    degenerate = point_with(
        {"kind": "isotropic"},
        target={"kind": "disc", "centre_m": [0, 0, 0], "normal": [0, 0, 1], "radius_m": 2.0},
        position_m=[0.0, 0.0, 1e-170],
    )
    with pytest.raises(ValueError, match="2\\^53"):
        photon_count(degenerate, seed=1)
    with pytest.raises(ValueError, match="2\\^53"):
        next(iter_chunks(degenerate, seed=1))


def test_complex_pulse_widths_must_be_finite():
    complex_pulse = {"kind": "complex", "flat_ns": 1.0, "gauss_rms_ns": 0.5, "decay_ns": 0.5}
    for name in ("flat_ns", "gauss_rms_ns", "decay_ns"):
        with pytest.raises(ValueError):
            emit_chunk(dict(FLASHER, pulse=dict(complex_pulse, **{name: math.inf})), seed=1, rows=5)
    assert expected_photons(dict(FLASHER, pulse=complex_pulse)) > 0.0
