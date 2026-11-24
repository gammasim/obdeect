"""Downward laser chain (sources -> propagation single scattering) against a numerical integral.

Chain: one point emission at height H above the observation level, axis straight down, cone 0,
355 nm, written with ``obdeect.sources.run.write``; then ``obdeect.propagation.propagate`` with
single scattering (weighted mode, 128 points per photon) and one telescope sphere at the origin.

Reference: the single-scattering integral of docs/propagation/scattering.md, integrated here in
altitude (independent code; the scattering table is parsed in this file). Per altitude step at
height z above the observation level, with P = (d, 0, z) and C the sphere centre:

    dW = N0 exp(-D(z)) [dtau_mol p_R(theta) + ssa dtau_aer p_HG(theta)] Omega T(P -> C)

with D the optical depth from the emission to P, theta the angle between the beam and P -> C,
Omega = 2 pi (1 - sqrt(1 - (R/D_PC)^2)) and T = exp(-od(P) D_PC / z) (od from the transmission
table, observation level to P).
"""

from pathlib import Path

import numpy as np
import pytest
from obdeect.sources import run

propagation = pytest.importorskip("obdeect.propagation")
stream = pytest.importorskip("obdeect.stream")

DATA = Path(__file__).resolve().parents[3] / "propagation" / "cpp" / "tests" / "data"
SCATTERING = DATA / "scattering_exponential.dat"
TRANSMISSION = DATA / "atm_trans_2000_1_10_0_0_2000.dat"
SITE_M = 2000.0
WAVELENGTH_NM = 355
EMISSION_HEIGHT_M = 100_000.0  # above the observation level
N0 = 1e15  # physical photons emitted
ROWS = 2000  # primary rows; bunch_size = N0 / ROWS
RADIUS_M = 12.5
N_POINTS = 128
SEED = 1
C_M_PER_NS = 0.299792458
DISCRETISATION_REL = 1.5e-5  # docs/propagation/scattering.md: relative error at 128 points
N_STEPS = 20000  # altitude steps of the reference (5 m)
DISTANCES_M = (100.0, 1000.0, 5000.0)  # horizontal offset of the beam from the telescope
# Scattering altitude bins [m a.s.l.] for the per-bin comparison.
BIN_EDGES_M = (2000.0, 3000.0, 6000.0, 15000.0, 40000.0, SITE_M + EMISSION_HEIGHT_M)


def scattering_table(path, wavelength_nm):
    """(altitude_km, tau_mol, tau_aer, ssa, g) of one wavelength of the scattering table."""
    rows = []
    for line in Path(path).read_text().splitlines():
        fields = line.split("#", 1)[0].split()
        if len(fields) == 6 and int(float(fields[0])) == wavelength_nm:
            rows.append([float(value) for value in fields[1:]])
    table = np.array(rows)
    assert len(table) > 1, "wavelength not in the scattering table"
    return table.T


def at_altitude(nodes_km, values, altitude_m):
    """Linear in log10(altitude), clamped to the edge nodes (as the scattering table)."""
    return np.interp(np.log10(altitude_m / 1000.0), np.log10(nodes_km), values)


def reference(d_m, transmission):
    """Altitude steps of the single-scattering integral: (z above observation level, weight)."""
    h_km, tau_mol, tau_aer, ssa_nodes, g_nodes = scattering_table(SCATTERING, WAVELENGTH_NM)
    z_edges = np.linspace(0.0, EMISSION_HEIGHT_M, N_STEPS + 1)
    h_edges = SITE_M + z_edges
    mol_edges = at_altitude(h_km, tau_mol, h_edges)
    aer_edges = at_altitude(h_km, tau_aer, h_edges)
    d_mol = mol_edges[:-1] - mol_edges[1:]  # vertical depth of each step (tau is depth above h)
    d_aer = aer_edges[:-1] - aer_edges[1:]

    z = 0.5 * (z_edges[:-1] + z_edges[1:])
    h = SITE_M + z
    top = np.array([SITE_M + EMISSION_HEIGHT_M])
    tau_top = at_altitude(h_km, tau_mol, top)[0] + at_altitude(h_km, tau_aer, top)[0]
    depth = at_altitude(h_km, tau_mol, h) + at_altitude(h_km, tau_aer, h) - tau_top
    ssa = at_altitude(h_km, ssa_nodes, h)
    g = at_altitude(h_km, g_nodes, h)

    dist = np.hypot(d_m, z)
    cos_theta = z / dist
    p_r = 3.0 / (16.0 * np.pi) * (1.0 + cos_theta**2)
    p_hg = (1.0 - g**2) / (4.0 * np.pi * (1.0 + g**2 - 2.0 * g * cos_theta) ** 1.5)
    omega = 2.0 * np.pi * (1.0 - np.sqrt(1.0 - (RADIUS_M / dist) ** 2))
    # The transmission table ends at 100 km; above 60 km the scattering depths do not change.
    od = transmission.optical_depth(np.full_like(h, WAVELENGTH_NM), np.minimum(h, 100_000.0))
    transmitted = np.exp(-od * dist / z)
    weight = np.exp(-depth) * (d_mol * p_r + ssa * d_aer * p_hg) * omega * transmitted
    return z, weight


def telescope_photons(path, telescope_id):
    """Concatenated photon columns of one telescope block, in stream order."""
    blocks = [
        record["columns"]
        for record in stream.Reader(str(path))
        if record["type"] == "photons" and record["telescope_id"] == telescope_id
    ]
    columns = {name: np.concatenate([block[name] for block in blocks]) for name in blocks[0]}
    columns["cz"] = np.concatenate([
        block.get("cz", -np.sqrt(np.maximum(0.0, 1.0 - block["cx"] ** 2 - block["cy"] ** 2)))
        for block in blocks
    ])
    return columns


def chain(directory, d_m):
    """Run sources and propagation; returns the telescope-1 photon columns."""
    description = {
        "run": {"number": 1, "seed": SEED, "observation_level_m": SITE_M, "frame": "sky"},
        "telescopes": [{"id": 1, "position_m": [0.0, 0.0, 0.0], "radius_m": RADIUS_M}],
        "events": [
            {
                "repeat": 1,
                "emissions": [
                    {
                        "kind": "point",
                        "photons": N0,
                        "bunch_size": N0 / ROWS,
                        "position_m": [d_m, 0.0, EMISSION_HEIGHT_M],
                        "axis": [0.0, 0.0, -1.0],
                        "angular": {"kind": "cone", "angle_deg": 0.0},
                        "spectrum": {"kind": "fixed", "wavelength_nm": float(WAVELENGTH_NM)},
                        "target": "none",
                    }
                ],
            }
        ],
    }
    summary = run.write(description, directory / "lidar_in.pst")
    assert summary["rows"] == ROWS
    propagation.propagate(
        directory / "lidar_in.pst",
        directory / "lidar_out.pst",
        transmission=TRANSMISSION,
        site_altitude_m=SITE_M,
        pathlen_m=0.0,  # no in-telescope path: the lidar chain tests the atmosphere only
        seed=SEED,
        scattering=SCATTERING,
        scattering_points=N_POINTS,
        mode="weighted",
    )
    return telescope_photons(directory / "lidar_out.pst", 1)


@pytest.fixture(scope="module")
def transmission():
    return propagation.TransmissionTable(str(TRANSMISSION), site_altitude_m=SITE_M)


@pytest.fixture(scope="module", params=DISTANCES_M, ids=lambda d: f"d{int(d)}m")
def case(request, transmission, tmp_path_factory):
    d_m = request.param
    photons = chain(tmp_path_factory.mktemp(f"d{int(d_m)}"), d_m)
    primary = slice(0, ROWS)
    scattered = slice(ROWS, len(photons["w"]))
    # Primaries come first in the block with cx = cy = 0 at the beam offset; scattered rows follow.
    assert np.all(photons["cx"][primary] == 0.0) and np.all(photons["cy"][primary] == 0.0)
    assert np.allclose(photons["x"][primary], d_m, rtol=0.0, atol=1e-3)
    assert np.all((photons["cx"][scattered] != 0.0) | (photons["cy"][scattered] != 0.0))
    out = {name: values[scattered].astype(np.float64) for name, values in photons.items()}
    z_ref, weight_ref = reference(d_m, transmission)
    return {"d": d_m, "scattered": out, "z_ref": z_ref, "weight_ref": weight_ref}


def test_total_scattered_weight_matches_integral(case):
    w = case["scattered"]["w"]
    chain_total = w.sum()
    reference_total = N0 * case["weight_ref"].sum()
    sigma = np.sqrt(np.sum(w**2) + (DISCRETISATION_REL * reference_total) ** 2)
    print(
        f"\nd={case['d']:.0f} m: rows {len(w)}, chain {chain_total:.6e}, "
        f"reference {reference_total:.6e}, ratio {chain_total / reference_total:.6f}, "
        f"sigma {sigma:.3e} ({sigma / reference_total:.2e} rel)"
    )
    assert reference_total > 0.0
    assert abs(chain_total - reference_total) <= 4.0 * sigma


def test_scattered_weight_per_altitude_bin(case):
    zem = case["scattered"]["zem"]
    w = case["scattered"]["w"]
    h_ref = case["z_ref"] + SITE_M
    for low, high in zip(BIN_EDGES_M[:-1], BIN_EDGES_M[1:]):
        chain_mask = (zem >= low) & (zem < high)
        ref_mask = (h_ref >= low) & (h_ref < high)
        chain_bin = w[chain_mask].sum()
        ref_bin = N0 * case["weight_ref"][ref_mask].sum()
        sigma = np.sqrt(np.sum(w[chain_mask] ** 2) + (DISCRETISATION_REL * ref_bin) ** 2)
        print(
            f"\nd={case['d']:.0f} m bin [{low / 1000:.0f}, {high / 1000:.0f}) km a.s.l.: "
            f"chain {chain_bin:.6e}, reference {ref_bin:.6e}, sigma {sigma:.3e}"
        )
        assert abs(chain_bin - ref_bin) <= 4.0 * sigma


def test_arrival_time_is_emission_to_point_plus_point_to_ground(case):
    photons = case["scattered"]
    d_m = case["d"]
    z_point = photons["zem"] - SITE_M
    # Ground position = output x, y (telescope at the origin); the point is at (d, 0, z).
    to_ground = np.sqrt((photons["x"] - d_m) ** 2 + photons["y"] ** 2 + z_point**2)
    expected = ((EMISSION_HEIGHT_M - z_point) + to_ground) / C_M_PER_NS
    delta = photons["t"] - expected
    # The ground point lies on the scattered ray: direction cosines towards it (cx, cy).
    assert np.allclose(photons["cx"], (photons["x"] - d_m) / to_ground, rtol=0.0, atol=1e-5)
    assert np.allclose(photons["cy"], photons["y"] / to_ground, rtol=0.0, atol=1e-5)
    print(
        f"\nd={case['d']:.0f} m: arrival-time error max |dt| {np.max(np.abs(delta)):.4f} ns, "
        f"mean {np.mean(delta):.4f} ns"
    )
    # Float32 altitude rounding is amplified by the slant path for horizontal rays.
    altitude_error = np.spacing(photons["zem"].astype(np.float32)).astype(np.float64) / 2
    time_error = altitude_error * (1 + 1 / np.abs(photons["cz"])) / C_M_PER_NS
    assert np.all(np.abs(delta) <= 0.1 + time_error)


def closest_approach(photons):
    """Distance to the telescope centre and travel time from the z=0 plane."""
    cx, cy = photons["cx"], photons["cy"]
    cz = photons["cz"]
    norm = np.sqrt(cx**2 + cy**2 + cz**2)
    cx, cy, cz = cx / norm, cy / norm, cz / norm
    along = -(photons["x"] * cx + photons["y"] * cy)
    distance = np.sqrt(
        (photons["x"] + along * cx) ** 2 + (photons["y"] + along * cy) ** 2 + (along * cz) ** 2
    )
    return distance, along / C_M_PER_NS


def test_arrival_spread_around_emission_height_over_c(case):
    photons = case["scattered"]
    w = photons["w"]
    _, time_shift = closest_approach(photons)
    delay = photons["t"] + time_shift - EMISSION_HEIGHT_M / C_M_PER_NS
    order = np.argsort(delay)
    cumulative = np.cumsum(w[order]) / w.sum()
    median, p95 = (delay[order][np.searchsorted(cumulative, q)] for q in (0.5, 0.95))
    print(f"\nd={case['d']:.0f} m: weighted delay median {median:.2f} ns, 95% {p95:.2f} ns")
    assert delay.min() >= -RADIUS_M / C_M_PER_NS - 0.1
    assert delay.max() <= (case["d"] + RADIUS_M) / C_M_PER_NS + 0.1


def test_scattered_rays_pass_through_sphere(case):
    distance, _ = closest_approach(case["scattered"])
    assert np.all(distance <= RADIUS_M + 1e-3)
