"""Single scattering: table lookup, zero optical depth, chunk independence, analytic yield."""

import math
from pathlib import Path

import numpy as np
import pytest
from dlufd import ScatteringTable, TransmissionTable, scatter

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
    return TransmissionTable(str(TRANSMISSION_PATH), site_altitude_m=OBSERVATION_LEVEL_M)


def _primaries(n, z_m=48000.0, weight=1.0, cx=None, cy=None, first=0):
    """Downward photons at height z_m above the observation level."""
    rng = np.random.default_rng(3)
    if cx is None:
        cx = np.zeros(n)
        cy = np.zeros(n)
    return {
        "x": np.zeros(n, dtype=np.float32),
        "y": np.zeros(n, dtype=np.float32),
        "z": np.full(n, z_m, dtype=np.float32),
        "cx": np.asarray(cx, dtype=np.float32),
        "cy": np.asarray(cy, dtype=np.float32),
        "t": rng.uniform(0.0, 1.0, n).astype(np.float32),
        "w": np.full(n, weight, dtype=np.float32),
        "wavelength": np.full(n, LAMBDA_NM, dtype=np.float32),
    }


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
