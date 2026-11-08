"""Print the install locations for CMake consumers: python -m gstrickd --cmake-dir."""

import argparse
import sys
from pathlib import Path

import gstrickd


def _locate(subdirectory: str, marker: str) -> Path | None:
    # Editable installs list the CMake install tree as a second package path.
    for root in gstrickd.__path__:
        candidate = Path(root) / subdirectory
        if (candidate / marker).is_file():
            return candidate
    return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m gstrickd",
        description="Print the directory of gstrickd's CMake package or C++ headers.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--cmake-dir",
        action="store_true",
        help="directory of gstrickd-config.cmake (use as -Dgstrickd_DIR=...)",
    )
    group.add_argument(
        "--include-dir",
        action="store_true",
        help="directory that contains the gstrickd/ headers (use as -I...)",
    )
    args = parser.parse_args(argv)

    if args.cmake_dir:
        path = _locate("cmake", "gstrickd-config.cmake")
    else:
        path = _locate("include", "gstrickd/stream.hpp")
    if path is None:
        print("gstrickd: no installed CMake package or C++ headers found", file=sys.stderr)
        return 1
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
