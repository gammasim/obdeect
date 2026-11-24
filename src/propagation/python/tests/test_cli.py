"""Command line of ``obdeect propagation run``: with extinction --pathlen is required, as it is for
the Python API (propagate); it has no default."""

import subprocess
import sys

from test_propagate import SITE_M, TABLE_PATH


def run_cli(tmp_path, *options):
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "obdeect.stream.cli",
            "propagation",
            "run",
            "--transmission",
            str(TABLE_PATH),
            "--site-altitude",
            str(SITE_M),
            "--seed",
            "1",
            *options,
            "-i",
            str(tmp_path / "in.pst"),
            "-o",
            str(tmp_path / "out.pst"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


def test_run_requires_pathlen_with_extinction(tmp_path):
    result = run_cli(tmp_path)
    assert result.returncode == 2
    assert "--pathlen is required with extinction" in result.stderr
    assert not (tmp_path / "out.pst").exists()


def test_run_rejects_pathlen_with_no_extinction(tmp_path):
    result = run_cli(tmp_path, "--no-extinction", "--pathlen", "0")
    assert result.returncode == 2
    assert "--no-extinction cannot be combined with --transmission, --pathlen" in result.stderr
