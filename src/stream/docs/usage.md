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
$ ps-cat --version                       # ps-cat <package version> (photon-stream format 1.0)
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
find_package(gstrickd 0.1 CONFIG REQUIRED)   # SameMinorVersion: any 0.1.x
target_link_libraries(<target> PRIVATE gstrickd::gstrickd_core)
```

Configure with `-Dgstrickd_DIR=$(python -m gstrickd --cmake-dir)`, or with
`-Dgstrickd_DIR=<prefix>/lib/cmake/gstrickd` for a `cmake --install` prefix.

The version is defined once. Wheels and sdists take it from the git tag (setuptools-scm,
passed by scikit-build-core as `SKBUILD_PROJECT_VERSION`). Other builds take the newest
`v*` tag, or 0.0.0 without one. CMake derives `project()` and the C++ `version.hpp` from it.

## Using gstrickd from other languages

`libgstrickd_c` exposes the reader and writer through a plain C ABI
(`gstrickd/gstrickd.h`, `gstrickd::gstrickd_c` in the CMake package), so Rust (`bindgen`
or hand-written `extern "C"`), Julia (`ccall`) or Go (cgo) can use it. Every returned
pointer, including the photon columns, belongs to the reader and is valid until the next
`gstrickd_reader_next()`; nothing is copied. Errors are negative codes with a message in
`gstrickd_last_error()`. Build scripts find the files with `python -m gstrickd --include-dir`,
`--lib-dir` and `--cmake-dir`.

A ray tracer in Python reads ground-frame photons and writes camera-frame photons:

```python
import numpy as np
from gstrickd.stage import transform

def trace(telescope_id, columns, metadata):
    # columns: x, y, cx, cy, t, w, ... as float32 arrays (ground frame, metres)
    x, y = columns["x"], columns["y"]
    hit = x**2 + y**2 < 12.0**2                  # your optics here
    out = {name: values[hit] for name, values in columns.items()}
    out["x"], out["y"] = out["x"] * 0.01, out["y"] * 0.01   # focal plane
    return out                                   # float32 arrays of equal length

transform("ground.ps", "camera.ps", trace, stage_id=0x444C0010, frame_out="camera",
          ledger={"mirror": {"state": "applied", "source": "my-model-v1"}})
```

`transform` copies all other records unchanged, sets the frame of the output chunks and
telescope records, and appends one `ledger` record per run (it refuses a factor that is
already `applied`).

The same in C (error handling shortened):

```c
#include <gstrickd/gstrickd.h>

gstrickd_reader *in;  gstrickd_writer *out;
gstrickd_reader_open("ground.ps", &in);
gstrickd_writer_open("camera.ps", "my-tracer", &out);
int type;
while (gstrickd_reader_next(in, &type) == GSTRICKD_OK) {
  if (type == GSTRICKD_RECORD_PHOTONS) {
    uint32_t telescope, rows;  uint64_t first;
    gstrickd_reader_photons(in, &telescope, NULL, &first, &rows, NULL);
    const float *columns[GSTRICKD_N_COLUMNS];
    for (int c = 0; c < GSTRICKD_N_COLUMNS; ++c)
      columns[c] = gstrickd_reader_column(in, c);        /* NULL if absent; trace here */
    gstrickd_writer_write_photons(out, telescope, GSTRICKD_FRAME_CAMERA, first, rows, columns);
  }                       /* copy other records with gstrickd_reader_meta_* / write_record */
}                         /* the loop ends with GSTRICKD_END, anything else is an error */
gstrickd_writer_close(out);  gstrickd_reader_close(in);
```

Link with `-L$(python -m gstrickd --lib-dir) -lgstrickd_c` or, in CMake,
`target_link_libraries(tracer PRIVATE gstrickd::gstrickd_c)`.

### Conformance streams

`cpp/tests/data/conformance/` holds `valid.ps` (every record type, plus a record of unknown
type that must be skipped), four invalid streams (bad magic, unknown major version,
truncated chunk, more than 2^16 rows) and `expected.json` with the decoded content and the
expected error. `make_conformance.py` regenerates them with the standard library only. An
independent implementation of [](format.md) should decode `valid.ps` to `expected.json` and
reject the other files.
