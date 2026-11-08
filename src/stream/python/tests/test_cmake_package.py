"""The CMake package: python -m gstrickd paths, and a consumer built against them.

Set GSTRICKD_CMAKE_DIR to test another install, e.g. a `cmake --install` prefix:
GSTRICKD_CMAKE_DIR=build/prefix/lib/cmake/gstrickd pytest python/tests/test_cmake_package.py
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
CONSUMER = ROOT / "cpp" / "tests" / "consumer"


def _module_output(option):
    result = subprocess.run(
        [sys.executable, "-m", "gstrickd", option],
        check=True,
        capture_output=True,
        text=True,
    )
    return Path(result.stdout.strip())


def test_cmake_dir_and_include_dir():
    assert (_module_output("--cmake-dir") / "gstrickd-config.cmake").is_file()
    assert (_module_output("--include-dir") / "gstrickd" / "stream.hpp").is_file()


def test_consumer_builds_and_runs(tmp_path):
    cmake = shutil.which("cmake")
    if cmake is None:
        pytest.skip("cmake not found")
    cmake_dir = os.environ.get("GSTRICKD_CMAKE_DIR") or str(_module_output("--cmake-dir"))
    build = tmp_path / "build"
    subprocess.run(
        [cmake, "-S", str(CONSUMER), "-B", str(build), f"-Dgstrickd_DIR={cmake_dir}"],
        check=True,
        capture_output=True,
    )
    subprocess.run([cmake, "--build", str(build)], check=True, capture_output=True)
    executable = build / ("consumer.exe" if os.name == "nt" else "consumer")
    subprocess.run([str(executable)], check=True)
