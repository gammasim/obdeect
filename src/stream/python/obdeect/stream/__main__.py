"""Print the install locations for CMake consumers: python -m obdeect.stream --cmake-dir."""

import argparse
import sys
from pathlib import Path

import obdeect.stream as obdeect_stream


def _locate(subdirectory: str, marker: str) -> Path | None:
    # Editable installs list the CMake install tree as a second package path.
    # marker is a file name or a glob pattern.
    for root in obdeect_stream.__path__:
        candidate = Path(root) / subdirectory
        if candidate.is_dir() and next(candidate.glob(marker), None) is not None:
            return candidate
    return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m obdeect.stream",
        description="Print the directory of the obdeect stream CMake package, headers or C library.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--cmake-dir",
        action="store_true",
        help="directory of obdeect_stream-config.cmake (use as -Dobdeect_stream_DIR=...)",
    )
    group.add_argument(
        "--include-dir",
        action="store_true",
        help="directory that contains the obdeect/stream/ headers, C and C++ (use as -I...)",
    )
    group.add_argument(
        "--lib-dir",
        action="store_true",
        help="directory of the C library libobdeect_stream_c (use as -L... and an rpath)",
    )
    args = parser.parse_args(argv)

    if args.cmake_dir:
        path = _locate("cmake", "obdeect_stream-config.cmake")
    elif args.lib_dir:
        path = _locate("lib", "*obdeect_stream_c.*")
    else:
        path = _locate("include", "obdeect/stream/stream.hpp")
    if path is None:
        print(
            "obdeect stream: no installed CMake package, headers or library found", file=sys.stderr
        )
        return 1
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
