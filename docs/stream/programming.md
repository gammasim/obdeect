# Programming interfaces

The format is available in C++, Python and C. Record types and fields are in [](format.md).

## C++

```cpp
#include "obdeect/stream/stream.hpp"

#include <fstream>

std::ofstream file("run.pst", std::ios::binary);
obdeect::stream::Writer writer(file, "my-module");           // writes the stream record
obdeect::stream::PhotonChunk chunk;
chunk.telescope_id = 1;
chunk.frame = obdeect::stream::Frame::telescope;
chunk.resize(3, obdeect::stream::required_columns);          // x y cx cy t w
chunk[obdeect::stream::Column::x] = {0.0F, 1.0F, 2.0F};      // fill the columns, then:
writer.write_chunk(chunk);
writer.finish();                                      // writes stream_end

std::ifstream in("run.pst", std::ios::binary);
obdeect::stream::Reader reader(in);
obdeect::stream::Record record;
while (reader.next(record)) {
  if (record.type == obdeect::stream::RecordType::photons) {
    // record.chunk.rows(), record.chunk[obdeect::stream::Column::w], ...
  }
}
```

## Python

```python
import numpy as np
import obdeect.stream

with obdeect.stream.Writer("run.pst", "my-module") as writer:
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

for record in obdeect.stream.Reader("run.pst"):
    print(record["type"])                             # stream, run, photons, blob, ...
    if record["type"] == "photons":
        w = record["columns"]["w"]                    # numpy float32 array
```

Metadata values: `int`, `float`, `str`, or a list or array of floats. Photon columns must be
float32 arrays. Names: `obdeect.stream.COLUMNS`, `obdeect.stream.FRAMES`. `Reader("-")` reads
stdin.

The module also provides `obdeect.stream.corsika_to_stream`, `stream_to_simtel` and
`stream_to_grisu`, the Python counterparts of the [commands](command-line.md).

## CMake package

`cmake --install` puts the headers and a CMake package in the install prefix. A wheel carries the
same files under `obdeect/stream/include` and `obdeect/stream/cmake`:

```console
$ python -m obdeect.stream --cmake-dir      # directory of obdeect_stream-config.cmake
$ python -m obdeect.stream --include-dir    # directory of the obdeect/stream/ headers
```

A consumer links the header-only core:

```cmake
find_package(obdeect stream 0.1 CONFIG REQUIRED)   # SameMinorVersion: any 0.1.x
target_link_libraries(<target> PRIVATE obdeect::stream)
```

Configure with `-Dobdeect_stream_DIR=$(python -m obdeect.stream --cmake-dir)`, or with
`-Dobdeect_stream_DIR=<prefix>/lib/cmake/obdeect_stream` for a `cmake --install` prefix. Other
consumers and build options: [](../user-guide/installation.md).

The version comes from the git tag (setuptools-scm, passed by scikit-build-core as
`SKBUILD_PROJECT_VERSION`). Other builds take the newest `v*` tag, or 0.0.0 without one.

## Other languages and custom stages

`libobdeect_stream_c` exposes the reader and writer through a plain C ABI
(`obdeect/stream/obdeect_stream.h`, `obdeect::stream_c` in the CMake package). Rust (`bindgen` or
`extern "C"`), Julia (`ccall`) and Go (cgo) can use it. Returned pointers, including photon
columns, belong to the reader and are valid until the next `obdeect_stream_reader_next()`; nothing
is copied. Errors are negative codes with a message in `obdeect_stream_last_error()`.

The optics step, which turns ground-frame photons into camera-frame photons, is
`obdeect ray-tracing trace-stream` ([](../ray-tracing/photon-stream.md)). A Python ray tracer
uses the same stage interface:

```python
import numpy as np
from obdeect.stream.stage import transform

def trace(telescope_id, columns, metadata):
    # columns: x, y, cx, cy, t, w, ... as float32 arrays (ground frame, metres)
    x, y = columns["x"], columns["y"]
    hit = x**2 + y**2 < 12.0**2                  # your optics here
    out = {name: values[hit] for name, values in columns.items()}
    out["x"], out["y"] = out["x"] * 0.01, out["y"] * 0.01   # focal plane
    return out                                   # float32 arrays of equal length

transform("ground.pst", "camera.pst", trace, stage_id=0x444C0010, frame_out="camera",
          efficiency_status={"mirror": {"state": "applied", "source": "my-model-v1"}})
```

`transform` copies all other records unchanged (records of unknown type are dropped), sets the
frame of the output chunks and telescope records, and appends one `efficiency_status` record per
run. It refuses a factor that is `applied` in any input efficiency status.

The same in C (error handling shortened):

```c
#include <obdeect/stream/obdeect_stream.h>

obdeect_stream_reader *in;  obdeect_stream_writer *out;
obdeect_stream_reader_open("ground.pst", &in);
obdeect_stream_writer_open("camera.pst", "my-tracer", &out);
int type;
while (obdeect_stream_reader_next(in, &type) == OBDEECT_STREAM_OK) {
  if (type == OBDEECT_STREAM_RECORD_PHOTONS) {
    uint32_t telescope, rows;  uint64_t first;
    obdeect_stream_reader_photons(in, &telescope, NULL, &first, &rows, NULL);
    const float *columns[OBDEECT_STREAM_N_COLUMNS];
    for (int c = 0; c < OBDEECT_STREAM_N_COLUMNS; ++c)
      columns[c] = obdeect_stream_reader_column(in, c);        /* NULL if absent; trace here */
    obdeect_stream_writer_write_photons(out, telescope, OBDEECT_STREAM_FRAME_CAMERA, first, rows, columns);
  }                       /* copy other records with obdeect_stream_reader_meta_* / write_record */
}                         /* the loop ends with OBDEECT_STREAM_END, anything else is an error */
obdeect_stream_writer_close(out);  obdeect_stream_reader_close(in);
```

Link with `-L$(python -m obdeect.stream --lib-dir) -lobdeect_stream_c`, or in CMake with
`target_link_libraries(tracer PRIVATE obdeect::stream_c)`.

## Conformance streams

`cpp/tests/data/conformance/` holds `valid.pst` (every record type, plus a record of unknown type
that must be skipped), four invalid streams (bad magic, unknown major version, truncated chunk,
more than 2^16 rows) and `expected.json` with the decoded content and the expected error.
`make_conformance.py` regenerates them with the standard library only. An independent
implementation of [](format.md) should decode `valid.pst` to `expected.json` and reject the others.

## End-to-end check

`validation/stream/validate_container.sh OUTDIR` is a manual end-to-end check of the EventIO codec. In the
`ghcr.io/gammasim/simtools-dev` container it runs CORSIKA, then
`obdeect stream from-corsika | obdeect stream to-simtel`, and compares the CORSIKA EventIO file
with the round-trip file. It then runs sim_telarray on both files with the same configuration and
seed and compares the outputs. The container, CORSIKA and sim_telarray are not part of the
automated tests.
