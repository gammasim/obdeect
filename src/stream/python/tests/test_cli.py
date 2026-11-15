"""The obdeect command: component groups, native programs and versions."""

import subprocess
import sys

import obdeect.stream as obdeect_stream
import pytest
from obdeect.stream import cli

OBDEECT = [sys.executable, "-m", "obdeect.stream.cli"]


def test_stream_group_is_registered():
    group = cli.groups()["stream"]
    assert set(group.commands) == {"cat", "tee", "from-corsika", "to-simtel", "to-grisu"}


@pytest.mark.parametrize("command", ["cat", "tee", "from-corsika", "to-simtel", "to-grisu"])
def test_native_programs_are_installed(command):
    path = cli.program_path("obdeect.stream", f"obdeect-stream-{command}")
    result = subprocess.run([path, "--version"], capture_output=True, text=True, check=True)
    assert result.stdout.startswith(f"obdeect-stream-{command} ") or result.stdout.startswith(
        f"obdeect stream {command} "
    )


def test_command_runs_the_native_program():
    result = subprocess.run(
        [*OBDEECT, "stream", "cat", "--version"], capture_output=True, text=True, check=True
    )
    assert result.stdout.startswith(f"obdeect stream cat {obdeect_stream.__version__[:3]}")


def test_version_lists_components_and_format(capsys):
    assert cli.main(["--version"]) == 0
    out = capsys.readouterr().out
    assert f"obdeect-stream {obdeect_stream.__version__}" in out
    assert f"photon-stream format {obdeect_stream.format_major}." in out


def test_missing_component_names_its_distribution(monkeypatch, capsys):
    monkeypatch.setattr(cli, "groups", lambda component=None: {})
    assert cli.main(["sources", "led"]) == 2
    assert "pip install obdeect-sources" in capsys.readouterr().err


def test_groups_can_load_one_component():
    assert set(cli.groups("stream")) == {"stream"}
    assert cli.groups("no-such-component") == {}


def test_command_loads_only_its_component(monkeypatch):
    # A command imports its own component only: the other components (ray-tracing imports NumPy)
    # are not loaded, which keeps every `obdeect <component> <command>` start-up short.
    calls = []
    monkeypatch.setattr(
        cli, "run_command", lambda group, name, command, argv: calls.append(name) or 0
    )
    monkeypatch.delitem(sys.modules, "obdeect.ray_tracing.cli", raising=False)
    assert cli.main(["stream", "cat", "file.pst"]) == 0
    assert calls == ["stream cat"]
    assert "obdeect.ray_tracing.cli" not in sys.modules


def test_unknown_command_is_an_error():
    result = subprocess.run([*OBDEECT, "stream", "nonsense"], capture_output=True, text=True)
    assert result.returncode == 2
    assert "invalid choice" in result.stderr
