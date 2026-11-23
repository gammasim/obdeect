"""Runs: write a run description as a photon-stream file, and check it without generating photons.

``write`` stores the normalized description (JSON, sorted keys) as the run metadata key
``run_description``, so the file says how it was made. ``check`` returns the candidates and
expected weights per emission (see ``_core.check_run``).

.. code-block:: python

    summary = run.write("run.toml", "run.pst", seed=3)  # events, rows, weight, expected
    numbers = run.check("run.toml")
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from obdeect.sources import _core, config


def _normalized(description: Mapping[str, Any] | str | Path) -> dict[str, Any]:
    """Normalized run from a description or a .toml / .json file name."""
    if isinstance(description, (str, Path)):
        description = config.load(description)
    return config.normalize(description)


def write(
    description: Mapping[str, Any] | str | Path,
    output: str | Path,
    *,
    seed: int | None = None,
    max_rows: float | None = None,
) -> dict[str, Any]:
    """Write a run to a photon-stream file (``-`` = stdout); returns the summary.

    ``seed`` overrides ``run.seed``. ``max_rows`` rejects a run whose expected candidate rows
    (all telescopes) exceed it, before anything is written. The summary holds ``events``,
    ``rows`` (all telescopes), ``weight`` (sum of weights written), ``expected`` (expected sum
    of weights of all emitted photons) and ``not_on_plane`` (photons of ``to_plane_z0`` emissions
    that move away from the plane z = 0 and are dropped; they are in ``expected`` but not in
    ``rows`` or ``weight``).
    """
    run = _normalized(description)
    if seed is not None:
        if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
            raise ValueError("seed must be an int in 0..2^64-1")
        run["run"]["seed"] = seed
    if max_rows is not None:
        if not math.isfinite(max_rows) or max_rows < 0:
            raise ValueError("max_rows must be finite and non-negative")
        rows = _core.check_run(run)["expected_candidates"]
        if rows > max_rows:
            raise ValueError(
                f"the run would write about {rows:.3g} photon rows (limit {max_rows:.3g}, "
                f"~{rows * 40 / 1e9:.3g} GB); raise bunch_size or the limit (--max-rows)"
            )
    run["run"]["metadata"]["run_description"] = json.dumps(run, sort_keys=True)
    return _core.write_run(run, str(output))


def check(description: Mapping[str, Any] | str | Path) -> dict[str, Any]:
    """Validated run with the candidates and expected weights of every emission.

    Keys: ``events``, ``telescopes``, ``emissions`` (one entry per emission: ``event``,
    ``emission``, ``repeat``, ``seed``, ``candidates``, ``expected_weight`` per occurrence),
    ``expected_candidates`` and ``expected_weight`` (run totals, all telescopes).
    """
    return _core.check_run(_normalized(description))
