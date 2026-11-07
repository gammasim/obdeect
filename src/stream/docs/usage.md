# Usage

A stream is written and read record by record. See [](format.md) for the record
types and fields.

## C++

```cpp
#include "gstrickd/stream.hpp"

#include <fstream>

std::ofstream file("run.ps", std::ios::binary);
gstrickd::Writer writer(file, "my-module");           // writes the stream record
gstrickd::PhotonChunk chunk;
chunk.telescope_id = 1;
chunk.frame = gstrickd::Frame::telescope;
chunk.resize(3, gstrickd::required_columns);          // x y cx cy t w
chunk[gstrickd::Column::x] = {0.0F, 1.0F, 2.0F};      // fill the columns, then:
writer.write_chunk(chunk);
writer.finish();                                      // writes stream_end

std::ifstream in("run.ps", std::ios::binary);
gstrickd::Reader reader(in);
gstrickd::Record record;
while (reader.next(record)) {
  if (record.type == gstrickd::RecordType::photons) {
    // record.chunk.rows(), record.chunk[gstrickd::Column::w], ...
  }
}
```

## Python

```python
import numpy as np
import gstrickd

with gstrickd.Writer("run.ps", "my-module") as writer:
    writer.write("run", {"run_number": 1, "telescope_x_m": [0.0, 10.0]})
    writer.write_photons(
        1,
        "telescope",                                  # or 2
        0,
        {"x": np.zeros(3, np.float32), "y": np.zeros(3, np.float32),
         "cx": np.zeros(3, np.float32), "cy": np.zeros(3, np.float32),
         "t": np.zeros(3, np.float32), "w": np.ones(3, np.float32)},
    )
    writer.write_blob(1204, 0, b"...")                # raw EventIO block

for record in gstrickd.Reader("run.ps"):
    print(record["type"])                             # stream, run, photons, blob, ...
    if record["type"] == "photons":
        w = record["columns"]["w"]                    # numpy float32 array
```

Metadata values: `int`, `float`, `str`, or a list or array of floats. Photon columns
must be float32 arrays; the names are `gstrickd.COLUMNS`, the frames are `gstrickd.FRAMES`.

`Reader("-")` reads from stdin.

## Command line

```console
$ ps-cat run1.ps run2.ps > all.ps        # one stream: one stream and one stream_end record
$ ps-cat --summary all.ps                # records per type, photons and weight per telescope
$ ps-cat run.ps | ps-tee copy.ps > out.ps  # copy.ps and out.ps get the same records
$ corsika2ps run.iact > run.ps           # CORSIKA IACT EventIO to a photon stream
$ ps2simtel --camera-map cam.txt run.ps > sim.simtel   # stream to sim_telarray input
$ corsika2ps < run.iact | ps2simtel > sim.simtel       # pipe; "-" reads stdin
$ ps-cat --version                       # ps-cat 0.0.1 (photon-stream format 1.0)
```

The camera map lists the camera-frame telescopes, one per line:

```text
# telescope_id mirror_class(1 single | 2 dual) camera_rotation_deg
2 1 0.0
```

Telescopes not listed are ground-frame (sim_telarray `Bypass_Optics=0`). Listed telescopes are
converted to the bypass convention and need `Bypass_Optics=1` (class 1) or `2` (class 2) in
sim_telarray; `ps2simtel` prints the value for each.

## CMake package

`cmake --install` puts the headers and a CMake package in the install prefix; a wheel
carries the same files under `gstrickd/include` and `gstrickd/cmake`:

```console
$ python -m gstrickd --cmake-dir      # directory of gstrickd-config.cmake
$ python -m gstrickd --include-dir    # directory of the gstrickd/ headers
```

A consumer (e.g. dlufd or dlamb) links the header-only core:

```cmake
find_package(gstrickd 0.0 CONFIG REQUIRED)   # SameMinorVersion
target_link_libraries(<target> PRIVATE gstrickd::gstrickd_core)
```

Configure with `-Dgstrickd_DIR=$(python -m gstrickd --cmake-dir)`, or with
`-Dgstrickd_DIR=<prefix>/lib/cmake/gstrickd` for a `cmake --install` prefix.

The version is defined once. Wheels and sdists take it from the git tag (setuptools-scm,
passed by scikit-build-core as `SKBUILD_PROJECT_VERSION`). Other builds take the newest
`v*` tag, or 0.0.0 without one. CMake derives `project()` and the C++ `version.hpp` from it.
