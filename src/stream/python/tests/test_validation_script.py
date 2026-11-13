"""Failure propagation in the container validation's required commands."""

import os
import shlex
import subprocess
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    ("corsika_status", "adapter_status", "expected"), [(42, 0, 42), (0, 43, 43)]
)
def test_required_commands_propagate_failure(tmp_path, corsika_status, adapter_status, expected):
    script = (Path(__file__).parents[4] / "validation/stream/validate_container.sh").read_text()
    fragment = script.split("start=$(date +%s)", 1)[1].split("cmp -l run.iact", 1)[0]
    corsika = tmp_path / "corsika"
    corsika.write_text(f"#!/bin/sh\nprintf partial > run.iact\nexit {corsika_status}\n")
    source = tmp_path / "obdeect-stream-from-corsika"
    source.write_text(f"#!/bin/sh\ncat\nexit {adapter_status}\n")
    sink = tmp_path / "obdeect-stream-to-simtel"
    sink.write_text("#!/bin/sh\ncat\n")
    for path in (corsika, source, sink):
        path.chmod(0o755)
    fragment = fragment.replace(
        "/workdir/simulation_software/corsika7/corsika_qgs2_urqmd_flat", shlex.quote(str(corsika))
    ).replace("/out/bin/obdeect-stream-from-corsika", shlex.quote(str(source)))
    fragment = fragment.replace("/out/bin/obdeect-stream-to-simtel", shlex.quote(str(sink)))
    (tmp_path / "card.txt").write_text("")
    result = subprocess.run(
        ["bash", "-c", "set -uo pipefail\nstart=0\n" + fragment],
        cwd=tmp_path,
        env={**os.environ, "LC_ALL": "C"},
        capture_output=True,
        check=False,
    )
    assert result.returncode == expected
