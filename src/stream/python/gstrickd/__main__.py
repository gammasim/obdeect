"""Print the install locations for CMake consumers: python -m gstrickd --cmake-dir."""

import argparse
import sys
from pathlib import Path

import gstrickd


def _locate(subdirectory: str, marker: str) -> Path | None:
    # Editable installs list the CMake install tree as a second package path.
    # marker is a file name or a glob pattern.
    for root in gstrickd.__path__:
        candidate = Path(root) / subdirectory
        if candidate.is_dir() and next(candidate.glob(marker), None) is not None:
            return candidate
    return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m gstrickd",
        description="Print the directory of gstrickd's CMake package, headers or C library.",
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
        help="directory that contains the gstrickd/ headers, C and C++ (use as -I...)",
    )
    group.add_argument(
        "--lib-dir",
        action="store_true",
        help="directory of the C library libgstrickd_c (use as -L... and an rpath)",
    )
    args = parser.parse_args(argv)

    if args.cmake_dir:
        path = _locate("cmake", "gstrickd-config.cmake")
    elif args.lib_dir:
        path = _locate("lib", "*gstrickd_c.*")
    else:
        path = _locate("include", "gstrickd/stream.hpp")
    if path is None:
        print("gstrickd: no installed CMake package, headers or library found", file=sys.stderr)
        return 1
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
