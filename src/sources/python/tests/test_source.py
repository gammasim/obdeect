"""Light sources from dictionaries: frames, distributions, chunks and validation."""

import math

import numpy as np
import pytest
from dlamb.source import (
    CHUNK_ROWS,
    COLUMNS,
    build_source,
    emit,
    emit_chunk,
    expected_photons,
    iter_chunks,
    photon_count,
)

FLASHER = {
    "frame": "camera",
    "photons": 4e5,
    "emitter": {
        "position_m": [0.0, 0.0, 1.5],
        "pointing": [0.0, 0.0, -1.0],
        "angular": {"kind": "filled_cone", "angle_deg": 11.0},
    },
    "spectrum": {"kind": "fixed", "wavelength_nm": 385.0},
    "pulse": {"kind": "tophat", "width_ns": 2.5},
    "target": {
        "kind": "disc",
        "centre_m": [0.0, 0.0, 0.0],
        "normal": [0.0, 0.0, 1.0],
        "radius_m": 1.2,
    },
}

STAR = {
    "frame": "telescope",
    "photons": 1e5,
    "star": {
        "towards_source": [0.01, 0.0, 1.0],
        "distance_m": math.inf,
        "disk_radius_m": 5.0,
        "disk_z_m": 0.0,
        "start_z_m": 30.0,
    },
}


def emitter_with(angular, photons=2e5, **extra):
    spec = {
        "frame": "sky",
        "photons": photons,
        "emitter": {
            "position_m": [0.0, 0.0, 10.0],
            "pointing": [0.0, 0.0, -1.0],
            "angular": angular,
        },
    }
    spec.update(extra)
    return spec


def cosine_to_axis(chunk, axis):
    axis = np.asarray(axis, dtype=float)
    return (chunk["cx"] * axis[0] + chunk["cy"] * axis[1] + chunk["cz"] * axis[2]).astype(float)


def test_columns_are_float32_and_directions_unit():
    chunk = emit_chunk(FLASHER, seed=3, rows=5000)
    assert tuple(chunk) == COLUMNS
    assert all(array.dtype == np.float32 for array in chunk.values())
    length = np.sqrt(chunk["cx"] ** 2 + chunk["cy"] ** 2 + chunk["cz"] ** 2)
    assert np.allclose(length, 1.0, atol=1e-6)
    assert np.all((chunk["w"] >= 0.0) & (chunk["w"] <= 1.0))
    # Positions are the emitter (x, y, z) and the target points lie on the disc plane.
    assert np.all(chunk["z"] == np.float32(1.5))
    assert np.all(chunk["wavelength"] == np.float32(385.0))


def test_disc_target_weights_sum_to_photons_reaching_the_disc():
    # The 11 deg cone (pointing at the disc) lands entirely on the 1.2 m disc: the weights
    # sum to the emitted number of photons (within the Poisson-free rounding of the count).
    photons = emit(FLASHER, seed=3)
    assert photons["w"].sum() == pytest.approx(4e5, rel=1e-3)
    assert len(photons["w"]) == photon_count(FLASHER, seed=3)


def test_star_plane_wave_starts_on_the_start_plane():
    # obdeect PointSource (infinite distance): rays start at z = start_z and are parallel
    # to -towards_source; a quarter of the disc (radius 2.5 of 5 m) is hit by 25 % of them.
    photons = emit(STAR, seed=2)
    assert np.allclose(photons["z"], 30.0)
    towards = np.array([0.01, 0.0, 1.0]) / np.hypot(0.01, 1.0)
    assert np.allclose(photons["cx"], -towards[0], atol=1e-6)
    assert np.allclose(photons["cz"], -towards[2], atol=1e-6)
    s = 30.0 / photons["cz"].astype(float)  # Ray parameter to the disc plane z = 0.
    hit = np.hypot(photons["x"] - s * photons["cx"], photons["y"] - s * photons["cy"])
    assert np.all(hit <= 5.0 + 1e-3)
    inside = np.count_nonzero(hit < 2.5) / len(hit)
    assert inside == pytest.approx(0.25, abs=4.0 * math.sqrt(0.25 * 0.75 / len(hit)))


def test_point_source_at_finite_distance_emits_from_its_position():
    spec = {
        "frame": "telescope",
        "photons": 20.0,
        "star": {
            "towards_source": [0.6, 0.0, 0.8],
            "distance_m": 20.0,
            "disk_radius_m": 1.0,
            "disk_z_m": 1.5,
            "start_z_m": 30.0,
        },
    }
    photons = emit(spec, seed=1)
    assert np.allclose(photons["x"], 12.0, atol=1e-4)
    assert np.allclose(photons["z"], 17.5, atol=1e-4)
    # t = (|start - source| - distance) / c with start = source: -distance / c.
    assert np.allclose(photons["t"], -20.0 / 0.299792458, atol=1e-3)


@pytest.mark.parametrize(
    ("angular", "expected_cosine"),
    [
        ({"kind": "cone", "angle_deg": 10.0}, math.cos(math.radians(10.0))),
        ({"kind": "lambertian"}, 2.0 / 3.0),  # E[cos theta] for intensity cos theta.
        ({"kind": "isotropic"}, 0.0),
    ],
)
def test_angular_profiles_mean_cosine(angular, expected_cosine):
    photons = emit(emitter_with(angular), seed=5)
    cosine = cosine_to_axis(photons, [0.0, 0.0, -1.0])
    sigma = max(np.std(cosine), 1e-12) / math.sqrt(len(cosine))
    assert cosine.mean() == pytest.approx(expected_cosine, abs=4.0 * sigma + 1e-6)


def test_gauss_pulse_and_exponential_moments():
    gauss = emit(
        emitter_with(
            {"kind": "isotropic"}, pulse={"kind": "gauss", "sigma_ns": 2.0, "mean_ns": 5.0}
        ),
        seed=4,
    )
    assert gauss["t"].mean() == pytest.approx(5.0, abs=4.0 * 2.0 / math.sqrt(len(gauss["t"])))
    assert gauss["t"].std() == pytest.approx(2.0, rel=0.02)
    decay = emit(
        emitter_with(
            {"kind": "isotropic"}, pulse={"kind": "exponential", "decay_ns": 3.0, "mean_ns": 1.0}
        ),
        seed=4,
    )
    assert decay["t"].mean() == pytest.approx(1.0, abs=4.0 * 3.0 / math.sqrt(len(decay["t"])))


def test_tabulated_spectrum_and_pulse_follow_their_tables():
    spec = emitter_with(
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
    photons = emit(emitter_with({"kind": "gauss", "sigma_deg": 5.0}, photons=1e5), seed=2)
    angle = np.degrees(np.arccos(np.clip(cosine_to_axis(photons, [0.0, 0.0, -1.0]), -1, 1)))
    # Two-dimensional r.m.s. of the Gaussian angular profile: sigma (LightEmission "Gauss").
    assert np.sqrt(np.mean(angle**2) / 2.0) == pytest.approx(5.0, rel=0.02)


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
    spec = dict(STAR, photons=2.5)
    counts = [photon_count(spec, seed=seed) for seed in range(4000)]
    assert set(counts) <= {2, 3}
    assert np.mean(counts) == pytest.approx(2.5, abs=4.0 * 0.5 / math.sqrt(len(counts)))
    # Candidates per photon of the 11 deg cone on the disc (obdeect Emitter formula):
    # I_max / (2 pi integral(I du)) * pi R^2 / h^2 with integral(I du) = 1 - cos(11 deg).
    per_photon = 1.0 / (2.0 * math.pi * (1.0 - math.cos(math.radians(11.0))))
    per_photon *= math.pi * 1.2**2 / 1.5**2
    assert expected_photons(FLASHER) == pytest.approx(4e5 * per_photon, rel=1e-12)


def test_sphere_target_photons_all_hit_the_sphere():
    spec = emitter_with(
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


def test_invalid_sources_are_rejected():
    def rejected(spec, match=None):
        with pytest.raises(ValueError, match=match):
            build_source(spec)

    rejected(dict(FLASHER, unknown=1), "unknown source key")
    rejected(dict(FLASHER, frame="ground"), "frame")
    rejected({k: v for k, v in FLASHER.items() if k != "photons"}, "photons")
    rejected(dict(FLASHER, photons=-1.0))
    rejected(dict(FLASHER, star=STAR["star"]), "exactly one")
    rejected({"frame": "sky", "photons": 1.0}, "exactly one")
    rejected(dict(FLASHER, spectrum={"kind": "nonsense"}), "unknown spectrum")
    rejected(dict(FLASHER, pulse={"kind": "gauss", "sigma_ns": 1.0, "delay": 2.0}))
    rejected(dict(FLASHER, pulse={"kind": "gauss", "mean_ns": 1.0}), "sigma_ns")
    rejected(dict(FLASHER, angular_extra=1), "unknown")
    rejected(
        emitter_with(
            {"kind": "cone", "angle_deg": 5.0},
            target={"kind": "disc", "centre_m": [0, 0, 0], "normal": [0, 0, 1], "radius_m": 1.0},
        ),
        "fixed-angle",
    )
    rejected(
        dict(
            STAR,
            target={"kind": "disc", "centre_m": [0, 0, 0], "normal": [0, 0, 1], "radius_m": 1.0},
        ),
        "only emitters",
    )
    rejected(dict(STAR, star=dict(STAR["star"], towards_source=[0.0, 0.0, 0.0])), "zero")
    with pytest.raises(ValueError, match="rows"):
        emit_chunk(FLASHER, seed=1, rows=CHUNK_ROWS + 1)
    with pytest.raises(ValueError, match="first_index"):
        emit_chunk(FLASHER, seed=1, first_index=-1)


def test_kind_is_allowed_only_in_profile_and_target_sections():
    def rejected(spec, match):
        with pytest.raises(ValueError, match=match):
            build_source(spec)

    rejected(dict(FLASHER, kind="disc"), "unknown source key")
    rejected(
        dict(FLASHER, emitter=dict(FLASHER["emitter"], kind="isotropic")), "unknown emitter key"
    )
    rejected(dict(STAR, star=dict(STAR["star"], kind="plane")), "unknown star key")
    # Accepted where a kind is expected: profiles and targets.
    assert build_source(FLASHER) is not None
    assert build_source(STAR) is not None


@pytest.mark.parametrize("angle", [-10.0, 190.0, math.nan, math.inf])
def test_cone_angles_are_degrees_in_0_180(angle):
    # Negative and > 180 angles used to alias to valid cones through 1 - cos.
    for kind in ("cone", "filled_cone"):
        with pytest.raises(ValueError, match="finite and in"):
            build_source(emitter_with({"kind": kind, "angle_deg": angle}))
    assert build_source(emitter_with({"kind": "cone", "angle_deg": 180.0})) is not None


@pytest.mark.parametrize("angles", [[0.0, -5.0, 10.0], [0.0, 200.0], [0.0, math.nan, 30.0]])
def test_angular_table_angles_are_checked_before_the_cosine(angles):
    table = {"kind": "table", "angle_deg": angles, "value": [1.0] * len(angles)}
    with pytest.raises(ValueError, match="finite and in"):
        build_source(emitter_with(table))
    ok = {"kind": "table", "angle_deg": [0.0, 90.0, 180.0], "value": [1.0, 0.5, 0.1]}
    assert build_source(emitter_with(ok)) is not None


def test_gauss_angular_width_must_be_representable():
    with pytest.raises(ValueError):
        build_source(emitter_with({"kind": "gauss", "sigma_deg": 1e-200}))


def test_photon_counts_must_be_finite_and_exact():
    # Counts beyond 2^53 are not exact in a double; iter_chunks would loop without end.
    with pytest.raises(ValueError):
        build_source(dict(FLASHER, photons=1e20))
    with pytest.raises(ValueError):
        build_source(dict(FLASHER, photons=math.inf))
    # A disc target at h = 1e-170 makes the candidate count overflow to inf.
    degenerate = emitter_with(
        {"kind": "isotropic"},
        target={"kind": "disc", "centre_m": [0, 0, 0], "normal": [0, 0, 1], "radius_m": 2.0},
    )
    degenerate["emitter"]["position_m"] = [0.0, 0.0, 1e-170]
    with pytest.raises(ValueError, match="2\\^53"):
        photon_count(degenerate, seed=1)
    with pytest.raises(ValueError, match="2\\^53"):
        next(iter_chunks(degenerate, seed=1))


def test_complex_pulse_widths_must_be_finite():
    complex_pulse = {"kind": "complex", "flat_ns": 1.0, "gauss_rms_ns": 0.5, "decay_ns": 0.5}
    for name, value in (("flat_ns", math.inf), ("gauss_rms_ns", math.inf), ("decay_ns", math.inf)):
        pulse = dict(complex_pulse, **{name: value})
        with pytest.raises(ValueError):
            build_source(dict(FLASHER, pulse=dict(pulse, delay_ns=0.0)))
    assert build_source(dict(FLASHER, pulse=dict(complex_pulse, delay_ns=0.0))) is not None
