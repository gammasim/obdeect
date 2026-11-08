"""Atmospheric extinction against sim_telarray reference values and photon statistics."""

import math
from pathlib import Path

import numpy as np
import pytest
from dlufd import TransmissionTable, apply_extinction

DATA = Path(__file__).resolve().parents[2] / "cpp" / "tests" / "data"
TABLE_PATH = DATA / "atm_trans_2000_1_10_0_0_2000.dat"
OBSERVATION_LEVEL_M = 2000.0
TELESCOPE_ALTITUDE_M = 2030.0
PATHLEN_M = 28.0

# (lambda_nm, zem_m, airmass, transmission) from the simtools Python port of sim_telarray's
# transmission lookup (simtools/src/simtools/simtel/shower_transmission.py), obs 2000 m,
# telescope 2030 m, pathlen 28 m; same values as cpp/tests/test_transmission.cpp.
REFERENCE = [
    (300, 2000.0, 1.0, 0.98806406601046559),
    (300, 2000.0, 1.7, 0.98375386766848683),
    (300, 2000.0, 2.5, 0.97885094506037973),
    (300, 2500.0, 1.0, 0.90783842908383527),
    (300, 2500.0, 1.7, 0.85185628219901666),
    (300, 2500.0, 2.5, 0.7920912361299689),
    (300, 3000.0, 1.0, 0.8336981623718448),
    (300, 3000.0, 1.7, 0.73699884847261166),
    (300, 3000.0, 2.5, 0.64014139543126103),
    (300, 4000.0, 1.0, 0.7216255829643391),
    (300, 4000.0, 1.7, 0.57660982796059346),
    (300, 4000.0, 2.5, 0.44620414981990048),
    (300, 8000.0, 1.0, 0.46867633888778365),
    (300, 8000.0, 1.7, 0.27684451467767029),
    (300, 8000.0, 2.5, 0.15168281275069989),
    (300, 20000.0, 1.0, 0.22152056318596161),
    (300, 20000.0, 1.7, 0.077438272082548013),
    (300, 20000.0, 2.5, 0.023296405593862524),
    (400, 2000.0, 1.0, 0.99479019142911018),
    (400, 2000.0, 1.7, 0.99289837070916542),
    (400, 2000.0, 2.5, 0.9907406948698555),
    (400, 2500.0, 1.0, 0.96035262500394203),
    (400, 2500.0, 1.7, 0.93517630561819554),
    (400, 2500.0, 2.5, 0.90721053607743829),
    (400, 3000.0, 1.0, 0.92949231080643413),
    (400, 3000.0, 1.7, 0.88466555857594431),
    (400, 3000.0, 2.5, 0.8360760239990529),
    (400, 4000.0, 1.0, 0.88484415612701106),
    (400, 4000.0, 1.7, 0.81364469141074114),
    (400, 4000.0, 2.5, 0.73926174530871436),
    (400, 8000.0, 1.0, 0.77522052303516964),
    (400, 8000.0, 1.7, 0.64980672574849396),
    (400, 8000.0, 2.5, 0.53112215759287196),
    (400, 20000.0, 1.0, 0.68153903988685849),
    (400, 20000.0, 1.7, 0.5220300606473498),
    (400, 20000.0, 2.5, 0.38490904004396648),
    (450, 2000.0, 1.0, 0.99585002816243451),
    (450, 2000.0, 1.7, 0.99433551018430177),
    (450, 2000.0, 2.5, 0.99260745268672734),
    (450, 2500.0, 1.0, 0.96885426807141428),
    (450, 2500.0, 1.7, 0.94894846579812375),
    (450, 2500.0, 2.5, 0.92669927585910872),
    (450, 3000.0, 1.0, 0.94505980048990335),
    (450, 3000.0, 1.7, 0.90967036119903621),
    (450, 3000.0, 2.5, 0.87084508928690763),
    (450, 4000.0, 1.0, 0.91165279142923306),
    (450, 4000.0, 1.7, 0.85568386680255726),
    (450, 4000.0, 2.5, 0.79591437198388915),
    (450, 8000.0, 1.0, 0.83071784993134223),
    (450, 8000.0, 1.7, 0.73059092533525782),
    (450, 8000.0, 2.5, 0.63085062504629652),
    (450, 20000.0, 1.0, 0.76435489719562222),
    (450, 20000.0, 1.7, 0.63416875458624955),
    (450, 20000.0, 2.5, 0.51230702687194263),
    (500, 2000.0, 1.0, 0.99664713490561485),
    (500, 2000.0, 1.7, 0.99542629946045724),
    (500, 2000.0, 2.5, 0.9940328900120935),
    (500, 2500.0, 1.0, 0.97466494081448352),
    (500, 2500.0, 1.7, 0.95839111327345694),
    (500, 2500.0, 2.5, 0.940124907990597),
    (500, 3000.0, 1.0, 0.95557452981351076),
    (500, 3000.0, 1.7, 0.92669848601867155),
    (500, 3000.0, 2.5, 0.8947642245184928),
    (500, 4000.0, 1.0, 0.92948322669951544),
    (500, 4000.0, 1.7, 0.88409582391534669),
    (500, 4000.0, 2.5, 0.83493196899756361),
    (500, 8000.0, 1.0, 0.86728683778110249),
    (500, 8000.0, 1.7, 0.78589653967824646),
    (500, 8000.0, 2.5, 0.70218910442835591),
    (500, 20000.0, 1.0, 0.81846755913236324),
    (500, 20000.0, 1.7, 0.71218232627378286),
    (500, 20000.0, 2.5, 0.60750642256193665),
    (600, 2000.0, 1.0, 0.9966093268583528),
    (600, 2000.0, 1.7, 0.99535324790559221),
    (600, 2000.0, 2.5, 0.9939196674784524),
    (600, 2500.0, 1.0, 0.97843574956029356),
    (600, 2500.0, 1.7, 0.96469440247482874),
    (600, 2500.0, 2.5, 0.94922615695124912),
    (600, 3000.0, 1.0, 0.96359889355951822),
    (600, 3000.0, 1.7, 0.93995816580564462),
    (600, 3000.0, 2.5, 0.91364955694379024),
    (600, 4000.0, 1.0, 0.94469499851512928),
    (600, 4000.0, 1.7, 0.90882564837441071),
    (600, 4000.0, 2.5, 0.86949673696933916),
    (600, 8000.0, 1.0, 0.90280368122902888),
    (600, 8000.0, 1.7, 0.84138249793051578),
    (600, 8000.0, 2.5, 0.77628683220750205),
    (600, 20000.0, 1.0, 0.87111336423902919),
    (600, 20000.0, 1.7, 0.79179321837256622),
    (600, 20000.0, 2.5, 0.70994647581709736),
    (700, 2000.0, 1.0, 0.98996269888848132),
    (700, 2000.0, 1.7, 0.98611900171347877),
    (700, 2000.0, 2.5, 0.98174447560465361),
    (700, 2500.0, 1.0, 0.95996055533928093),
    (700, 2500.0, 1.7, 0.93585394369890074),
    (700, 2500.0, 2.5, 0.90904390524019851),
    (700, 3000.0, 1.0, 0.94024003159666369),
    (700, 3000.0, 1.7, 0.90340641237797692),
    (700, 3000.0, 2.5, 0.8630743959829078),
    (700, 4000.0, 1.0, 0.91916825955625225),
    (700, 4000.0, 1.7, 0.86925830168738327),
    (700, 4000.0, 2.5, 0.8155281120288671),
    (700, 8000.0, 1.0, 0.88455809748627934),
    (700, 8000.0, 1.7, 0.81435193182397081),
    (700, 8000.0, 2.5, 0.74091313640882439),
    (700, 20000.0, 1.0, 0.86736527349324111),
    (700, 20000.0, 1.7, 0.78762731009385378),
    (700, 20000.0, 2.5, 0.70543412511710468),
]

ENVELOPE_WAVELENGTH = [200.0, 800.0]
ENVELOPE_VALUE = [0.9, 0.3]
WAVELENGTH_RANGE = (250.0, 700.0)


@pytest.fixture(scope="module")
def table():
    return TransmissionTable(str(TABLE_PATH), OBSERVATION_LEVEL_M)


def test_reference_transmission(table):
    od_tel, od_focus = table.telescope_depths(TELESCOPE_ALTITUDE_M, PATHLEN_M)
    points = np.array(REFERENCE, dtype=np.float64)
    lam = points[:, 0]
    index = lam.astype(int)
    values = table.transmission2(lam, points[:, 1], points[:, 2], od_tel[index], od_focus[index])
    relative = np.abs(values - points[:, 3]) / points[:, 3]
    assert relative.max() <= 1e-6


def test_reference_transmission_with_od2focus_from_table(table):
    # sim_telarray's od2focus for a 28 m path is tau(h + 28 m) - tau(h), taken from the table.
    points = np.array(REFERENCE, dtype=np.float64)
    lam, zem, airmass, expected = points.T
    od_tel = np.array([table.optical_depth(lm, TELESCOPE_ALTITUDE_M) for lm in lam])
    od_focus = (
        np.array([table.optical_depth(lm, TELESCOPE_ALTITUDE_M + PATHLEN_M) for lm in lam]) - od_tel
    )
    values = table.transmission2(lam, zem, airmass, od_tel, od_focus)
    relative = np.abs(values - expected) / expected
    assert relative.max() <= 1e-6


def test_constant_optical_depth(tmp_path):
    path = tmp_path / "constant.dat"
    path.write_text("# H2= 2.000, H1=    3.000     4.000\n400 0.5 0.5\n", encoding="utf-8")
    constant = TransmissionTable(str(path), 2000.0)
    for airmass in (1.0, 1.5, 2.5):
        result = constant.transmission(400.0, 3500.0, airmass)
        assert result == pytest.approx(math.exp(-0.5 * airmass), rel=1e-14)
    assert constant.optical_depth(400.0, 5000.0) == pytest.approx(0.5, rel=1e-14)


def test_table_rejects_wrong_altitude():
    with pytest.raises(ValueError):
        TransmissionTable(str(TABLE_PATH), 1835.0)


def _photons(n, rng):
    return {
        "x": rng.uniform(-50, 50, n).astype(np.float32),
        "y": rng.uniform(-50, 50, n).astype(np.float32),
        "cx": rng.uniform(-0.2, 0.2, n).astype(np.float32),
        "cy": rng.uniform(-0.2, 0.2, n).astype(np.float32),
        "t": np.zeros(n, dtype=np.float32),
        "zem": rng.uniform(2500.0, 9000.0, n).astype(np.float32),
        "w": rng.choice([1.0, 2.5, 3.0, 5.5], n).astype(np.float32),
    }


def _apply(table, columns, weighted, first_index=0):
    return apply_extinction(
        columns,
        table,
        TELESCOPE_ALTITUDE_M,
        ENVELOPE_WAVELENGTH,
        ENVELOPE_VALUE,
        seed=20261107,
        event=5,
        telescope=2,
        first_index=first_index,
        weighted=weighted,
        wavelength_range=WAVELENGTH_RANGE,
        pathlen_m=PATHLEN_M,
    )


def test_thinned_count_matches_weighted_mean(table):
    columns = _photons(200000, np.random.default_rng(1))
    weighted = _apply(table, columns, weighted=True)
    thinned = _apply(table, columns, weighted=False)
    p = weighted["p"].astype(np.float64)
    w = columns["w"].astype(np.float64)
    whole = np.floor(w)
    frac = w - whole
    mean = np.sum(w * p)
    variance = np.sum(whole * p * (1 - p) + frac * p * (1 - frac * p))
    assert abs(len(thinned["x"]) - mean) <= 4 * math.sqrt(variance)
    assert np.all(thinned["w"] == 1.0)
    assert weighted["w"].sum() == pytest.approx(mean, rel=1e-5)
    assert set(thinned) == set(weighted)


def test_wavelengths_assigned_in_range(table):
    columns = _photons(20000, np.random.default_rng(2))
    weighted = _apply(table, columns, weighted=True)
    wavelength = weighted["wavelength"]
    assert np.all(wavelength == np.floor(wavelength))
    assert wavelength.min() >= WAVELENGTH_RANGE[0]
    assert wavelength.max() <= WAVELENGTH_RANGE[1]


def test_thinned_output_keeps_columns_when_all_rejected(table):
    # Envelope 0 rejects every photon; the columns are still returned, with zero rows.
    columns = _photons(50, np.random.default_rng(4))
    kwargs = dict(
        table=table,
        telescope_altitude_m=TELESCOPE_ALTITUDE_M,
        envelope_wavelength=ENVELOPE_WAVELENGTH,
        envelope_value=[0.0, 0.0],
        seed=20261107,
        event=5,
        telescope=2,
        first_index=0,
        wavelength_range=WAVELENGTH_RANGE,
        pathlen_m=PATHLEN_M,
    )
    weighted = apply_extinction(columns, weighted=True, **kwargs)
    thinned = apply_extinction(columns, weighted=False, **kwargs)
    assert set(thinned) == set(weighted) == set(columns) | {"wavelength", "p"}
    for name, values in thinned.items():
        assert len(values) == 0, name
        assert values.dtype == np.float32, name
    assert np.all(weighted["p"] == 0.0)


def test_empty_input_keeps_columns(table):
    columns = _photons(0, np.random.default_rng(5))
    result = _apply(table, columns, weighted=False)
    assert set(result) == set(columns) | {"wavelength", "p"}
    assert all(len(values) == 0 for values in result.values())


@pytest.mark.parametrize("bad", [400.5, np.nan, np.inf, -1.0])
def test_non_integral_or_non_finite_wavelength_rejected(table, bad):
    columns = _photons(10, np.random.default_rng(6))
    columns["wavelength"] = np.full(10, 400.0, dtype=np.float32)
    columns["wavelength"][3] = bad
    with pytest.raises(ValueError):
        _apply(table, columns, weighted=False)


@pytest.mark.parametrize("lambda_nm", [1000.0, 1e10])
def test_wavelength_at_or_above_max_is_dropped(table, lambda_nm):
    # sim_telarray never traces photons at or above MAX_LAMBDA; no int conversion overflows.
    columns = _photons(100, np.random.default_rng(7))
    columns["wavelength"] = np.full(100, lambda_nm, dtype=np.float32)
    thinned = _apply(table, columns, weighted=False)
    assert len(thinned["x"]) == 0
    weighted = _apply(table, columns, weighted=True)
    assert np.all(weighted["p"] == 0.0)
    assert np.all(weighted["w"] == 0.0)


def test_non_finite_cz_rejected(table):
    columns = _photons(10, np.random.default_rng(8))
    columns["cz"] = -np.sqrt(1.0 - columns["cx"].astype(np.float64) ** 2 - columns["cy"] ** 2)
    columns["cz"] = columns["cz"].astype(np.float32)
    _apply(table, columns, weighted=False)  # consistent cz is accepted
    columns["cz"][4] = np.nan
    with pytest.raises(ValueError):
        _apply(table, columns, weighted=False)


def test_non_finite_wavelength_range_rejected(table):
    columns = _photons(10, np.random.default_rng(9))
    with pytest.raises(ValueError):
        apply_extinction(
            columns,
            table,
            TELESCOPE_ALTITUDE_M,
            ENVELOPE_WAVELENGTH,
            ENVELOPE_VALUE,
            seed=1,
            event=1,
            telescope=1,
            first_index=0,
            weighted=False,
            wavelength_range=(250.0, np.inf),
            pathlen_m=PATHLEN_M,
        )


def test_invalid_input_rejected(table):
    columns = _photons(10, np.random.default_rng(3))
    columns["cx"][0] = np.float32(0.9)
    columns["cy"][0] = np.float32(0.9)
    with pytest.raises(ValueError):
        _apply(table, columns, weighted=False)
    columns = _photons(10, np.random.default_rng(3))
    columns["unknown"] = columns["w"]
    with pytest.raises(ValueError):
        _apply(table, columns, weighted=False)


@pytest.mark.parametrize("weighted", [True, False])
def test_split_input_matches_whole_input_bit_for_bit(table, weighted):
    # 200000 rows exceed one photon chunk (2^16); the binding slices internally. Splitting the
    # input by hand with matching first_index must give identical output.
    n = 200000
    columns = _photons(n, np.random.default_rng(10))
    first_index = 1000
    whole = _apply(table, columns, weighted=weighted, first_index=first_index)
    bounds = [0, 12345, 99999, n]
    parts = []
    for start, stop in zip(bounds[:-1], bounds[1:]):
        piece = {name: values[start:stop] for name, values in columns.items()}
        parts.append(_apply(table, piece, weighted=weighted, first_index=first_index + start))
    assert set(whole) == set(parts[0])
    for name, values in whole.items():
        joined = np.concatenate([part[name] for part in parts])
        assert values.dtype == joined.dtype, name
        assert values.tobytes() == joined.tobytes(), name
