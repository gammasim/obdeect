"""End-to-end chain: light source -> atmosphere -> optics -> sim_telarray input.

    obdeect sources point-source | obdeect propagation run | obdeect ray-tracing trace-stream
        | obdeect stream to-simtel

Needs all four components installed (pip install obdeect, or ./src/* from a checkout).
Checks that the components work together through the photon stream, not physics accuracy
(each component has its own tests).
"""

import subprocess
import sys
from pathlib import Path

import obdeect.stream as obdeect_stream
import pytest

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src" / "ray_tracing" / "python" / "tests"))
from telescopes import segmented  # noqa: E402

TRANSMISSION = REPOSITORY / "src/propagation/cpp/tests/data/atm_trans_2000_1_10_0_0_2000.dat"
OBDEECT = [sys.executable, "-m", "obdeect.stream.cli"]
CAMERA_FRAME = 3

SOURCE = ["sources", "point-source", "-z", "100000", "-n", "5e10", "-b", "10", "-r", "300"]
SOURCE += ["-h", "2000", "--seed", "3"]
ATMOSPHERE = ["propagation", "run", "--transmission", str(TRANSMISSION), "--site-altitude"]
ATMOSPHERE += ["2000", "--wavelength-range", "300", "600", "--seed", "5"]
OPTICS = ["ray-tracing", "trace-stream", "--seed", "7", "--pointing-zenith-deg", "0"]
OPTICS += ["--pointing-azimuth-deg", "0"]


def obdeect(*args, stdin=None):
    return subprocess.run(
        [*OBDEECT, *map(str, args)], input=stdin, capture_output=True, check=True
    ).stdout


@pytest.fixture(scope="module")
def model(tmp_path_factory):
    from obdeect.ray_tracing import OpticalModel

    records, files = segmented()
    path = tmp_path_factory.mktemp("model") / "model.json"
    OpticalModel.from_records(records, files.__getitem__, configuration_seed=1).save(path)
    return path


@pytest.fixture(scope="module")
def chain(model):
    """The stream after each step."""
    source = obdeect(*SOURCE)
    atmosphere = obdeect(*ATMOSPHERE, "-i", "-", "-o", "-", stdin=source)
    camera = obdeect(*OPTICS, "--model", model, stdin=atmosphere)
    return {"source": source, "atmosphere": atmosphere, "camera": camera}


def records(data, tmp_path, name):
    path = tmp_path / f"{name}.pst"
    path.write_bytes(data)
    return list(obdeect_stream.Reader(str(path)))


def photon_weight(recs):
    return sum(float(r["columns"]["w"].sum()) for r in recs if r["type"] == "photons")


def test_each_step_writes_a_complete_stream(chain, tmp_path):
    for name, data in chain.items():
        types = [r["type"] for r in records(data, tmp_path, name)]
        assert types[0] == "stream" and types[-1] == "run_end", name
        assert "photons" in types, name


def test_photons_reach_the_camera_frame(chain, tmp_path):
    camera = records(chain["camera"], tmp_path, "camera")
    telescopes = [r for r in camera if r["type"] == "telescope"]
    assert telescopes and all(t["metadata"]["frame"] == CAMERA_FRAME for t in telescopes)
    source = photon_weight(records(chain["source"], tmp_path, "source"))
    atmosphere = photon_weight(records(chain["atmosphere"], tmp_path, "atmosphere"))
    detected = photon_weight(camera)
    # Extinction and the optics each remove light, but never all of it.
    assert 0 < detected < atmosphere < source


def test_efficiency_status_records_the_atmosphere(chain, tmp_path):
    status = [
        r["metadata"]
        for r in records(chain["camera"], tmp_path, "c")
        if r["type"] == "efficiency_status"
    ]
    assert any(s.get("atmosphere.state") == "applied" for s in status)


def test_chain_is_deterministic(chain, model):
    source = obdeect(*SOURCE)
    assert source == chain["source"]
    atmosphere = obdeect(*ATMOSPHERE, "-i", "-", "-o", "-", stdin=source)
    assert atmosphere == chain["atmosphere"]
    assert obdeect(*OPTICS, "--model", model, stdin=atmosphere) == chain["camera"]


def test_camera_photons_convert_to_sim_telarray_input(chain, tmp_path):
    camera_map = tmp_path / "cam.txt"
    camera_map.write_text("1 1 0.0\n")
    simtel = tmp_path / "chain.simtel"
    obdeect(
        "stream", "to-simtel", "--camera-map", camera_map, "-o", simtel, "-", stdin=chain["camera"]
    )
    assert simtel.stat().st_size > len(chain["camera"]) // 2
