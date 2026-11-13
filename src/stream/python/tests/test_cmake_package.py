"""The CMake package: python -m obdeect.stream paths, and a consumer built against them.

Set OBDEECT_STREAM_CMAKE_DIR to test another install, e.g. a `cmake --install` prefix:
OBDEECT_STREAM_CMAKE_DIR=build/prefix/lib/cmake/obdeect_stream pytest python/tests/test_cmake_package.py
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
        [sys.executable, "-m", "obdeect.stream", option],
        check=True,
        capture_output=True,
        text=True,
    )
    return Path(result.stdout.strip())


def test_cmake_dir_and_include_dir():
    assert (_module_output("--cmake-dir") / "obdeect_stream-config.cmake").is_file()
    assert (_module_output("--include-dir") / "obdeect" / "stream" / "stream.hpp").is_file()
    assert (_module_output("--include-dir") / "obdeect" / "stream" / "obdeect_stream.h").is_file()
    assert any(_module_output("--lib-dir").glob("*obdeect_stream_c.*"))


def test_consumer_builds_and_runs(tmp_path):
    cmake = shutil.which("cmake")
    if cmake is None:
        pytest.skip("cmake not found")
    cmake_dir = os.environ.get("OBDEECT_STREAM_CMAKE_DIR") or str(_module_output("--cmake-dir"))
    build = tmp_path / "build"
    subprocess.run(
        [cmake, "-S", str(CONSUMER), "-B", str(build), f"-Dobdeect_stream_DIR={cmake_dir}"],
        check=True,
        capture_output=True,
    )
    subprocess.run([cmake, "--build", str(build)], check=True, capture_output=True)
    for name in ("consumer", "consumer_c"):
        executable = build / (name + (".exe" if os.name == "nt" else ""))
        subprocess.run([str(executable)], check=True)
