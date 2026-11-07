"""Atmospheric profile, refraction and timing through the Python API."""

import math
from pathlib import Path

import dlufd
import numpy as np
import pytest

PROFILE = (
    Path(__file__).resolve().parents[2]
    / "cpp"
    / "tests"
    / "data"
    / ("atmprof_ecmwf_north_winter_fixed.dat")
)
OBSERVATION_LEVEL_M = 2180.0
SPEED_OF_LIGHT_M_PER_NS = 0.299792458


@pytest.fixture(scope="module")
def profile():
    return dlufd.AtmosphereProfile(PROFILE)


@pytest.fixture(scope="module")
def tables(profile):
    return dlufd.RefractionTables(profile, OBSERVATION_LEVEL_M)


def test_profile_reads_corsika_table(profile):
    assert profile.levels == 55
    assert profile.bottom_m == 0.0
    assert profile.top_m == 120000.0


def test_vectorised_lookups_match_scalar_and_keep_shape(profile):
    heights = np.linspace(0.0, 119000.0, 24).reshape(4, 6)
    density = profile.density(heights)
    assert density.shape == heights.shape
    for index in np.ndindex(heights.shape):
        assert density[index] == pytest.approx(profile.density(float(heights[index])), rel=0)
    assert np.all(np.diff(profile.thickness(np.linspace(0.0, 119000.0, 500))) < 0.0)
    assert float(profile.refractive_index_minus_one(OBSERVATION_LEVEL_M)) == pytest.approx(
        2.251919e-4, rel=1e-5
    )


def test_group_index_matches_obdeect(profile):
    # obdeect propagation.ambient_group_index() on the same table (see test_profile.cpp).
    assert profile.group_index(2180.0) == pytest.approx(1.000225191893913, abs=1e-12)
    assert profile.group_index(1000.0) == pytest.approx(1.000255944992766, abs=1e-12)
    assert profile.group_index(60000.0) == pytest.approx(1.0000000583180002, abs=1e-12)


def test_vertical_photon_is_not_displaced(tables):
    height = np.array([3000.0, 10000.0])
    zero = np.zeros(2)
    cx, cy, x, y, t, bent = dlufd.refract(tables, height, zero, zero, zero, zero, zero)
    assert np.all(bent)
    assert np.all(cx == 0.0) and np.all(cy == 0.0)
    assert np.all(x == 0.0) and np.all(y == 0.0)
    vacuum_ns = (height - OBSERVATION_LEVEL_M) / SPEED_OF_LIGHT_M_PER_NS
    air_ns = np.array([
        dlufd.air_delay_ns(dlufd.AtmosphereProfile(PROFILE), 2180.0, h, 0.0) for h in height
    ])
    np.testing.assert_allclose(t - vacuum_ns, air_ns, rtol=0, atol=1e-6)


def test_refraction_displaces_and_bends_towards_vertical(tables):
    zenith = np.radians([20.0, 40.0, 60.0])
    cx = np.sin(zenith)
    cy = np.zeros_like(cx)
    n = cx.size
    new_cx, new_cy, x, y, t, bent = dlufd.refract(
        tables, np.full(n, 20000.0), cx, cy, np.zeros(n), np.zeros(n), np.zeros(n)
    )
    assert bent.dtype == bool and np.all(bent)
    assert np.all(new_cx < cx) and np.all(new_cx > 0.0)
    assert np.all(x < 0.0) and np.all(np.diff(x) < 0.0)
    assert np.all(np.abs(x) < 50.0)
    assert np.all(t > 0.0)


def test_upward_direction_is_left_unchanged(tables):
    one = np.array([1.1])
    zero = np.zeros(1)
    cx, cy, x, y, t, bent = dlufd.refract(tables, np.array([20000.0]), one, zero, zero, zero, zero)
    assert not bent[0]
    assert cx[0] == 1.1 and x[0] == 0.0 and t[0] == 0.0


def test_refract_rejects_mismatched_inputs(tables):
    good = np.zeros(3)
    with pytest.raises(ValueError):
        dlufd.refract(tables, good, good, good, good, good, np.zeros(2))
    with pytest.raises(ValueError):
        dlufd.refract(tables, np.zeros((3, 1)), good, good, good, good, good)


def test_time_helpers():
    assert dlufd.propagation_time_ns(1000.0, 1.0002) == pytest.approx(
        1000.0 * 1.0002 / SPEED_OF_LIGHT_M_PER_NS, rel=1e-14
    )
    delay = dlufd.air_delay_ns(dlufd.AtmosphereProfile(PROFILE), 2180.0, 10000.0, 0.0)
    assert delay == pytest.approx(4.0219, abs=1e-4)
    assert math.isfinite(delay)


def test_missing_profile_is_an_error(tmp_path):
    with pytest.raises(ValueError):
        dlufd.AtmosphereProfile(tmp_path / "missing.dat")
