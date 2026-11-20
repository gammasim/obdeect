"""Command line of ``obdeect sources``: dispatch, --set parsing, file or preset, error format.

``obdeect.sources.run`` (the extension) is replaced by a fake module, and the presets by a fake
preset module, so these tests need neither the extension nor the real presets.
"""

import argparse
import json
import sys
import types

import pytest
from obdeect.sources import cli, presets

MINIMAL = {
    "run": {"observation_level_m": 2150.0, "seed": 1},
    "telescopes": [{"radius_m": 6.0}],
    "events": [
        {
            "repeat": 1,
            "emissions": [
                {
                    "kind": "point",
                    "photons": 10.0,
                    "position_m": [0, 0, 100],
                    "spectrum": {"kind": "fixed", "wavelength_nm": 400.0},
                }
            ],
        }
    ],
}


@pytest.fixture
def fake_run(monkeypatch):
    """Fake ``obdeect.sources.run`` recording its calls."""
    calls = {"write": [], "check": []}
    module = types.ModuleType("obdeect.sources.run")

    def write(config, output, *, max_rows=None):
        calls["write"].append({"config": config, "output": output, "max_rows": max_rows})
        return {"events": 2, "rows": 7, "weight": 1.5, "expected": 2.0}

    def check(config):
        calls["check"].append(config)
        return {
            "events": 2,
            "emissions": [{"event": 0, "emission": 0, "kind": "point", "expected": 1.5}],
        }

    module.write = write
    module.check = check
    monkeypatch.setitem(sys.modules, "obdeect.sources.run", module)
    return calls


@pytest.fixture
def fake_preset(monkeypatch):
    """A preset with one option of its own (``--size``, which also defines ``--seed``)."""
    module = types.ModuleType("obdeect.sources.presets.fake_mod")
    module.NAME = "fake"
    module.HELP = "fake preset for tests"
    seen = {}

    def add_arguments(parser):
        parser.add_argument("--size", type=float, required=True, metavar="M")
        parser.add_argument("--seed", type=int, default=1, metavar="N")

    def build(args):
        seen["args"] = args
        return {**MINIMAL, "run": {**MINIMAL["run"], "seed": args.seed}, "size": args.size}

    module.add_arguments = add_arguments
    module.build = build
    monkeypatch.setitem(presets.PRESETS, "fake", "fake_mod")
    monkeypatch.setitem(sys.modules, "obdeect.sources.presets.fake_mod", module)
    return seen


def write_json(tmp_path, data=MINIMAL, name="run.json"):
    path = tmp_path / name
    path.write_text(json.dumps(data))
    return str(path)


# --- --set parsing ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("text", "path", "value"),
    [
        ("run.seed=3", ["run", "seed"], 3),
        ("events[0].repeat=5", ["events", 0, "repeat"], 5),
        ("telescopes[0].radius_m=12", ["telescopes", 0, "radius_m"], 12),
        ("events[1].emissions[2].photons=1e6", ["events", 1, "emissions", 2, "photons"], 1e6),
        ("target=auto", ["target"], "auto"),
        ("flag=true", ["flag"], True),
        ("values=[1, 2.5]", ["values"], [1, 2.5]),
        ("note=a=b", ["note"], "a=b"),
        ("empty=", ["empty"], ""),
    ],
)
def test_parse_set_values_and_paths(text, path, value):
    assert cli.parse_set(text) == (path, value)


@pytest.mark.parametrize("text", ["run.seed", "=3", "a..b", "a[x]=1", "[]=1", "a[0"])
def test_parse_set_rejects_bad_input(text):
    with pytest.raises(ValueError, match="--set"):
        cli.parse_set(text)


def test_apply_set_changes_nested_values():
    description = json.loads(json.dumps(MINIMAL))
    cli.apply_set(description, ["events", 0, "repeat"], 5, "events[0].repeat=5")
    cli.apply_set(description, ["telescopes", 0, "radius_m"], 12, "telescopes[0].radius_m=12")
    cli.apply_set(description, ["run", "metadata", "note"], "x", "run.metadata.note=x")
    assert description["events"][0]["repeat"] == 5
    assert description["telescopes"][0]["radius_m"] == 12
    assert description["run"]["metadata"] == {"note": "x"}


def test_apply_set_reports_index_out_of_range():
    with pytest.raises(ValueError, match="index 3 out of range"):
        cli.apply_set(json.loads(json.dumps(MINIMAL)), ["events", 3, "repeat"], 1, "e[3].repeat=1")


def test_apply_set_rejects_index_on_a_table():
    with pytest.raises(ValueError, match="needs a list"):
        cli.apply_set({"run": {}}, ["run", 0], 1, "run[0]=1")


# --- file or preset detection ----------------------------------------------------------------


def test_file_source_does_not_use_presets(tmp_path, monkeypatch, fake_run):
    def no_preset(name):
        raise AssertionError("presets must not be used for a file")

    monkeypatch.setattr(presets, "preset", no_preset)
    path = write_json(tmp_path)
    assert cli.run([path, "-o", str(tmp_path / "out.pst")]) == 0
    assert fake_run["write"][0]["config"]["run"]["seed"] == 1


def test_preset_source_does_not_load_files(tmp_path, monkeypatch, fake_run, fake_preset):
    def no_load(path):
        raise AssertionError("a preset name is not a file")

    monkeypatch.setattr(cli.config, "load", no_load)
    assert cli.run(["fake", "--size", "3", "-o", str(tmp_path / "out.pst")]) == 0
    assert fake_preset["args"].size == 3


# --- dispatch --------------------------------------------------------------------------------


def test_run_file_writes_with_overrides(tmp_path, fake_run, capsys):
    path = write_json(tmp_path)
    out = str(tmp_path / "out.pst")
    code = cli.run([
        path,
        "-o",
        out,
        "--seed",
        "4",
        "--set",
        "run.seed=9",
        "--set",
        "events[0].repeat=5",
    ])
    assert code == 0
    call = fake_run["write"][0]
    assert call["output"] == out
    assert call["config"]["run"]["seed"] == 9  # --set after --seed
    assert call["config"]["events"][0]["repeat"] == 5
    assert "events 2, rows 7, sum of weights 1.5, expected 2" in capsys.readouterr().err


def test_seed_option_alone_sets_run_seed(tmp_path, fake_run):
    assert cli.run([write_json(tmp_path), "-o", "x.pst", "--seed", "4"]) == 0
    assert fake_run["write"][0]["config"]["run"]["seed"] == 4


def test_overrides_do_not_change_the_loaded_file(tmp_path, fake_run):
    path = write_json(tmp_path)
    cli.run([path, "-o", "x.pst", "--set", "run.seed=9"])
    assert json.loads(open(path).read())["run"] == MINIMAL["run"]


def test_run_without_output_writes_to_stdout(tmp_path, fake_run):
    assert cli.run([write_json(tmp_path)]) == 0
    assert fake_run["write"][0]["output"] == "-"
    assert fake_run["write"][0]["max_rows"] == 1e8


def test_run_dash_output_is_stdout(tmp_path, fake_run):
    assert cli.run([write_json(tmp_path), "-o", "-", "--max-rows", "5"]) == 0
    assert fake_run["write"][0]["output"] == "-"
    assert fake_run["write"][0]["max_rows"] == 5.0


def test_preset_options_reach_build(fake_run, fake_preset, tmp_path):
    assert cli.run(["fake", "--size", "2.5", "--seed", "8", "-o", str(tmp_path / "o.pst")]) == 0
    assert fake_preset["args"].size == 2.5
    assert fake_run["write"][0]["config"]["run"]["seed"] == 8


def test_check_human_table(tmp_path, fake_run, capsys):
    assert cli.check([write_json(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "events: 2" in out
    assert "kind" in out and "expected" in out and "1.5" in out
    assert fake_run["check"]


def test_check_json_prints_normalized_run(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(cli.config, "normalize", lambda c: {"normalized": c["run"]["seed"]})
    assert cli.check([write_json(tmp_path), "--set", "run.seed=3", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == {"normalized": 3}


def test_check_json_with_real_normalize(tmp_path, capsys):
    assert cli.check([write_json(tmp_path), "--json"]) == 0
    normalized = json.loads(capsys.readouterr().out)
    assert normalized["run"]["seed"] == 1.0
    assert normalized["events"][0]["repeat"] == 1.0


def test_check_does_not_accept_seed(tmp_path, fake_run, capsys):
    assert cli.check([write_json(tmp_path), "--seed", "3"]) == 1
    assert "unrecognized arguments: --seed" in capsys.readouterr().err
    assert not fake_run["check"]


# --- help ------------------------------------------------------------------------------------


def test_run_help_lists_presets_with_help(fake_preset, capsys):
    assert cli.run(["--help"]) == 0
    out = capsys.readouterr().out
    assert "fake" in out and "fake preset for tests" in out


def test_preset_help_shows_its_options(fake_preset, capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.run(["fake", "--help"])
    assert exit_info.value.code == 0
    assert "--size" in capsys.readouterr().out


# --- errors ----------------------------------------------------------------------------------


def _one_error_line(capsys, command):
    err = capsys.readouterr().err
    assert err.count("\n") == 1
    assert "Traceback" not in err
    assert err.startswith(f"obdeect sources {command}: ")
    return err


def test_unknown_preset_is_one_line(capsys, fake_run):
    assert cli.run(["nope"]) == 1
    err = _one_error_line(capsys, "run")
    assert "unknown preset 'nope'" in err


def test_missing_source_is_one_line(capsys):
    assert cli.check([]) == 1
    _one_error_line(capsys, "check")


def test_source_must_come_first(capsys):
    assert cli.run(["--set", "run.seed=1", "x.json"]) == 1
    _one_error_line(capsys, "run")


def test_missing_preset_option_is_one_line(capsys, fake_preset):
    assert cli.run(["fake", "-o", "x.pst"]) == 1
    err = _one_error_line(capsys, "run")
    assert "--size" in err


def test_unknown_option_is_one_line(tmp_path, capsys, fake_run):
    assert cli.run([write_json(tmp_path), "--bogus"]) == 1
    err = _one_error_line(capsys, "run")
    assert "unrecognized arguments: --bogus" in err


def test_bad_set_is_one_line(tmp_path, capsys, fake_run):
    assert cli.run([write_json(tmp_path), "--set", "run.seed"]) == 1
    _one_error_line(capsys, "run")


def test_missing_file_is_one_line(tmp_path, capsys, fake_run):
    assert cli.run([str(tmp_path / "missing.toml")]) == 1
    _one_error_line(capsys, "run")


def test_value_error_from_write_is_one_line(tmp_path, capsys, fake_run, monkeypatch):
    module = sys.modules["obdeect.sources.run"]

    def write(config, output, *, max_rows=None):
        raise ValueError("events[0].emissions[0].photons: must be >= 0")

    monkeypatch.setattr(module, "write", write)
    assert cli.run([write_json(tmp_path), "-o", "x.pst"]) == 1
    err = _one_error_line(capsys, "run")
    assert "photons: must be >= 0" in err


# --- registration ----------------------------------------------------------------------------


def test_group_registers_run_and_check_as_python_commands():
    assert set(cli.group.commands) == {"run", "check"}
    for command in cli.group.commands.values():
        assert command.function is not None
        assert command.program is None


def test_namespace_is_not_used_as_option_store():
    # Options are read with getattr defaults, so a bare namespace also works.
    args = argparse.Namespace()
    assert cli._overrides(json.loads(json.dumps(MINIMAL)), args)["run"]["seed"] == 1
