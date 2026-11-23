"""Single scattering: table lookup, zero optical depth, chunk independence, analytic yield."""

import math
from pathlib import Path

import numpy as np
import pytest
from obdeect.propagation import ScatteringTable, TransmissionTable, scatter

DATA = Path(__file__).resolve().parents[2] / "cpp" / "tests" / "data"
SCATTERING_PATH = DATA / "scattering_exponential.dat"
TRANSMISSION_PATH = DATA / "atm_trans_2000_1_10_0_0_2000.dat"
OBSERVATION_LEVEL_M = 2000.0
LAMBDA_NM = 355.0
TELESCOPE = (1, 3000.0, 0.0, 0.0, 12.5)  # id, x, y, z above the observation level, radius [m]


@pytest.fixture(scope="module")
def table():
    return ScatteringTable(str(SCATTERING_PATH))


@pytest.fixture(scope="module")
def transmission():
    return TransmissionTable(
        str(TRANSMISSION_PATH),
        site_altitude_m=OBSERVATION_LEVEL_M,
        wavelength_interpolation="floor",
    )


def _primaries(n, z_m=48000.0, weight=1.0, cx=None, cy=None, first=0, cz=None):
    """Downward photons at height z_m above the observation level (cz: optional signed column)."""
    rng = np.random.default_rng(3)
    if cx is None:
        cx = np.zeros(n)
        cy = np.zeros(n)
    columns = {
        "x": np.zeros(n, dtype=np.float32),
        "y": np.zeros(n, dtype=np.float32),
        "z": np.full(n, z_m, dtype=np.float32),
        "cx": np.asarray(cx, dtype=np.float32),
        "cy": np.asarray(cy, dtype=np.float32),
        "t": rng.uniform(0.0, 1.0, n).astype(np.float32),
        "w": np.full(n, weight, dtype=np.float32),
        "wavelength": np.full(n, LAMBDA_NM, dtype=np.float32),
    }
    if cz is not None:
        columns["cz"] = np.asarray(cz, dtype=np.float32)
    return columns


def test_table_profile_at_node(table):
    # Node at 2 km of the example: tau_mol = tau_R0 exp(-2/8), tau_aer = tau_A0 exp(-2/1.5).
    tau_mol, tau_aer, ssa, g = table.profile(np.array([LAMBDA_NM]), np.array([2000.0]))
    tau_r0 = 0.00865 * (LAMBDA_NM / 1000.0) ** -4.05
    tau_a0 = 0.1 * (LAMBDA_NM / 550.0) ** -1.3
    assert tau_mol[0] == pytest.approx(tau_r0 * math.exp(-2.0 / 8.0), rel=1e-6)
    assert tau_aer[0] == pytest.approx(tau_a0 * math.exp(-2.0 / 1.5), rel=1e-6)
    assert ssa[0] == pytest.approx(0.9)
    assert g[0] == pytest.approx(0.65)


def test_zero_optical_depth_gives_no_photons(tmp_path, transmission):
    lines = [f"{wl} {h} 0 0 1 0" for wl in (354, 355, 356) for h in (1.0, 10.0, 60.0)]
    path = tmp_path / "zero.dat"
    path.write_text("\n".join(lines) + "\n")
    zero = ScatteringTable(str(path))
    out = scatter(_primaries(10), zero, transmission, [TELESCOPE], 1, 1, 0, 8)
    assert all(len(values) == 0 for values in out[TELESCOPE[0]].values())


def test_output_structure_and_bounds(table, transmission):
    out = scatter(_primaries(20, weight=2.0), table, transmission, [TELESCOPE], 7, 1, 0, 8)
    photons = out[TELESCOPE[0]]
    assert set(photons) == {"x", "y", "z", "cx", "cy", "cz", "t", "zem", "w", "wavelength"}
    n = len(photons["w"])
    assert n > 0
    assert all(v.dtype == np.float32 and len(v) == n for v in photons.values())
    assert np.all(photons["w"] > 0.0) and np.all(photons["w"] <= 2.0)
    assert np.allclose(photons["cx"] ** 2 + photons["cy"] ** 2 + photons["cz"] ** 2, 1.0, atol=1e-5)
    assert np.all(photons["cz"] < 0.0)  # downward towards the telescope
    # zem is the absolute altitude of the scattering point: z plus the observation level.
    assert np.allclose(photons["zem"], photons["z"] + OBSERVATION_LEVEL_M, atol=0.01)
    assert np.all(photons["wavelength"] == np.float32(LAMBDA_NM))


def test_chunk_split_is_identical(table, transmission):
    cx = np.linspace(-0.2, 0.2, 90)
    cy = np.linspace(0.1, -0.1, 90)
    full = _primaries(90, cx=cx, cy=cy)
    whole = scatter(full, table, transmission, [TELESCOPE], 11, 2, 0, 16)[TELESCOPE[0]]
    parts = []
    for begin, end in ((0, 30), (30, 90)):
        piece = {k: v[begin:end] for k, v in full.items()}
        parts.append(scatter(piece, table, transmission, [TELESCOPE], 11, 2, begin, 16)[1])
    for name in whole:
        joined = np.concatenate([p[name] for p in parts])
        assert np.array_equal(joined, whole[name]), name


def test_yield_agrees_with_analytic_integral(table, transmission):
    """Sum of scattered weights per primary against the single-scattering integral (numpy)."""
    h = np.linspace(OBSERVATION_LEVEL_M + 1e-3, 48000.0 + OBSERVATION_LEVEL_M, 400001)
    dist = np.hypot(TELESCOPE[1], h - OBSERVATION_LEVEL_M)
    mu = (h - OBSERVATION_LEVEL_M) / dist
    tau_mol, tau_aer, ssa, g = table.profile(np.full_like(h, LAMBDA_NM), h)
    dh = 0.5
    up = table.profile(np.full_like(h, LAMBDA_NM), h + dh)
    down = table.profile(np.full_like(h, LAMBDA_NM), h - dh)
    beta_mol = down[0] - up[0]
    beta_aer = down[1] - up[1]
    rayleigh = 3.0 / (16.0 * np.pi) * (1.0 + mu**2)
    hg = (1 - g**2) / (4 * np.pi * (1 + g**2 - 2 * g * mu) ** 1.5)
    phase = beta_mol * rayleigh + ssa * beta_aer * hg
    ratio = TELESCOPE[4] / dist
    solid = 2 * np.pi * (1 - np.sqrt(1 - ratio**2))
    od = transmission.optical_depth(np.full_like(h, LAMBDA_NM), h) - transmission.optical_depth(
        np.array([LAMBDA_NM]), np.array([OBSERVATION_LEVEL_M])
    )
    transmitted = np.exp(-od * dist / (h - OBSERVATION_LEVEL_M))
    top = table.profile(np.array([LAMBDA_NM]), np.array([48000.0 + OBSERVATION_LEVEL_M]))
    attenuation = np.exp(-((tau_mol + tau_aer) - (top[0] + top[1])))
    integrand = phase * solid * transmitted * attenuation
    analytic = np.trapezoid(integrand, h)

    n = 1000
    totals = np.array(
        [
            scatter(_primaries(1), table, transmission, [TELESCOPE], 5, 3, i, 128)[TELESCOPE[0]][
                "w"
            ].sum()
            for i in range(n)
        ],
        dtype=np.float64,
    )
    stderr = totals.std(ddof=1) / math.sqrt(n)
    assert abs(totals.mean() - analytic) <= 4.0 * stderr
    assert analytic > 0.0


TELESCOPE_ABOVE_GROUND = (1, 3000.0, 0.0, 500.7, 12.5)  # centre 500.7 m above the observation level


def _lidar_integral(table, transmission, telescope, h0, upward, h_lo, h_hi):
    """Single-scattering integral of a vertical primary (weight 1) into one sphere (numpy).

    Upward primaries: from h0 up to h_hi, attenuated by tau(h0) - tau(h); downward: from h0 down,
    attenuated by tau(h) - tau(h0). The sphere centre is the aim point. T and Omega as in
    scattering.hpp.
    """
    h_c = OBSERVATION_LEVEL_M + telescope[3]
    h = np.linspace(h_lo, h_hi, 400001)
    assert np.all(np.abs(h - h_c) > 0.0)  # the integrand is 0/0 at the centre height
    dist = np.hypot(telescope[1], h - h_c)
    mu = (1.0 if upward else -1.0) * (h_c - h) / dist
    lam = np.full_like(h, LAMBDA_NM)
    tau_mol, tau_aer, ssa, g = table.profile(lam, h)
    up = table.profile(lam, h + 0.5)
    down = table.profile(lam, h - 0.5)
    beta_mol = down[0] - up[0]
    beta_aer = down[1] - up[1]
    rayleigh = 3.0 / (16.0 * np.pi) * (1.0 + mu**2)
    hg = (1 - g**2) / (4 * np.pi * (1 + g**2 - 2 * g * mu) ** 1.5)
    phase = beta_mol * rayleigh + ssa * beta_aer * hg
    ratio = telescope[4] / dist
    solid = 2 * np.pi * (1 - np.sqrt(1 - ratio**2))
    od = np.abs(
        transmission.optical_depth(lam, h) - transmission.optical_depth(lam, np.full_like(h, h_c))
    )
    transmitted = np.exp(-od * dist / np.abs(h - h_c))
    h0_profile = table.profile(np.array([LAMBDA_NM]), np.array([h0]))
    depth_profile = tau_mol + tau_aer
    depth_top = h0_profile[0][0] + h0_profile[1][0]
    attenuation = (
        np.exp(-(depth_top - depth_profile)) if upward else np.exp(-(depth_profile - depth_top))
    )
    return np.trapezoid(phase * solid * transmitted * attenuation, h)


def test_upward_primary_yield_agrees_with_lidar_integral(table, transmission):
    """Upward laser from 1 m above the observation level into a sphere 500.7 m above it.

    The scattered weight per primary against the single-scattering integral over the whole upward
    path (up to 60 km, the table's top). Some scattering points are below the sphere centre: their
    scattered photons travel upwards.
    """
    h0 = OBSERVATION_LEVEL_M + 1.0
    analytic = _lidar_integral(table, transmission, TELESCOPE_ABOVE_GROUND, h0, True, h0, 60000.0)

    n = 1000
    totals = np.zeros(n)
    below = 0
    upward_rows = 0
    for i in range(n):
        out = scatter(
            _primaries(1, z_m=1.0, cz=np.ones(1)),
            table,
            transmission,
            [TELESCOPE_ABOVE_GROUND],
            5,
            3,
            i,
            128,
        )[TELESCOPE_ABOVE_GROUND[0]]
        totals[i] = out["w"].sum()
        below += int(np.sum(out["zem"] < OBSERVATION_LEVEL_M + 500.7))
        upward_rows += int(np.sum((out["zem"] < OBSERVATION_LEVEL_M + 500.7) & (out["cz"] > 0.0)))
    stderr = totals.std(ddof=1) / math.sqrt(n)
    assert below > 0
    assert upward_rows > 0
    assert analytic > 0.0
    assert abs(totals.mean() - analytic) <= 4.0 * stderr


def test_time_is_at_scattering_point(table, transmission):
    """t = t_in + (emission -> point) / c; the path from the point to the telescope is not added."""
    primaries = _primaries(20, z_m=48000.0, weight=1.0)
    primaries["t"] = np.full(20, 5.0, dtype=np.float32)
    photons = scatter(primaries, table, transmission, [TELESCOPE], 9, 1, 0, 16)[TELESCOPE[0]]
    assert len(photons["t"]) > 0
    # Vertical beam from h = 50000 m a.s.l. (z 48000 m above the observation level).
    expected = (
        5.0 + (48000.0 + OBSERVATION_LEVEL_M - photons["zem"].astype(np.float64)) / 0.299792458
    )
    assert np.allclose(photons["t"], expected, rtol=0.0, atol=0.05)


def test_rejects_missing_columns_and_non_finite_values(table, transmission):
    for name in ("x", "y", "cx", "cy", "t", "w"):
        primaries = _primaries(5)
        del primaries[name]
        with pytest.raises(ValueError):
            scatter(primaries, table, transmission, [TELESCOPE], 1, 1, 0, 8)
    for name, value in (("x", np.nan), ("y", np.inf), ("t", np.nan)):
        primaries = _primaries(5)
        primaries[name][2] = value
        with pytest.raises(ValueError):
            scatter(primaries, table, transmission, [TELESCOPE], 1, 1, 0, 8)


def test_profile_rejects_bad_altitude_and_wavelength(table):
    with pytest.raises(ValueError):
        table.profile(np.array([LAMBDA_NM]), np.array([-100.0]))
    with pytest.raises(ValueError):
        table.profile(np.array([LAMBDA_NM]), np.array([0.0]))
    with pytest.raises(ValueError):
        table.profile(np.array([np.nan]), np.array([2000.0]))
    # Huge wavelengths are clamped to the table edges (no integer overflow).
    edge = table.profile(np.array([1e300, -1e300]), np.array([2000.0, 2000.0]))
    last = table.profile(np.array([float(table.last_wavelength_nm)]), np.array([2000.0]))
    first = table.profile(np.array([float(table.first_wavelength_nm)]), np.array([2000.0]))
    assert edge[0][0] == last[0][0]
    assert edge[0][1] == first[0][0]
