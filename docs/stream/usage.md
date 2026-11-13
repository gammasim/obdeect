# Usage

Record types and fields: [](format.md). Installation: [](../installation.md).

## C++

```cpp
#include "obdeect/stream/stream.hpp"

#include <fstream>

std::ofstream file("run.ps", std::ios::binary);
obdeect::stream::Writer writer(file, "my-module");           // writes the stream record
obdeect::stream::PhotonChunk chunk;
chunk.telescope_id = 1;
chunk.frame = obdeect::stream::Frame::telescope;
chunk.resize(3, obdeect::stream::required_columns);          // x y cx cy t w
chunk[obdeect::stream::Column::x] = {0.0F, 1.0F, 2.0F};      // fill the columns, then:
writer.write_chunk(chunk);
writer.finish();                                      // writes stream_end

std::ifstream in("run.ps", std::ios::binary);
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

with obdeect.stream.Writer("run.ps", "my-module") as writer:
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

for record in obdeect.stream.Reader("run.ps"):
    print(record["type"])                             # stream, run, photons, blob, ...
    if record["type"] == "photons":
        w = record["columns"]["w"]                    # numpy float32 array
```

Metadata values: `int`, `float`, `str`, or a list or array of floats. Photon columns must be
float32 arrays. Names: `obdeect.stream.COLUMNS`, `obdeect.stream.FRAMES`. `Reader("-")` reads
stdin.

## Command line

```console
$ obdeect stream cat run1.ps run2.ps > all.ps          # one stream, one stream and one stream_end record
$ obdeect stream cat --summary all.ps                  # records per type, photons and weight per telescope
$ obdeect stream cat run.ps | obdeect stream tee copy.ps > out.ps   # copy.ps and out.ps get the same records
$ obdeect stream from-corsika run.iact > run.ps        # CORSIKA IACT EventIO to a photon stream
$ obdeect stream to-simtel --camera-map cam.txt run.ps > sim.simtel  # sim_telarray input
$ obdeect stream to-simtel --atmprof atmprof1.dat light.ps > sim.simtel  # streams without CORSIKA input
$ obdeect stream from-corsika < run.iact | obdeect stream to-simtel > sim.simtel  # pipe; "-" is stdin
$ obdeect stream to-grisu --telescope-config TLLOC.dat run.ps > run.grisu  # GrISU photon text
$ obdeect stream cat --version                         # package version (photon-stream format 1.0)
```

### Camera map and telescope frames

The camera map lists the camera-frame telescopes, one per line:

```text
# telescope_id mirror_class(1 single | 2 dual) camera_rotation_deg
2 1 0.0
```

Telescopes not listed are ground-frame (sim_telarray `Bypass_Optics=0`). Listed telescopes with
camera-frame photons are converted to the bypass convention and need `Bypass_Optics=1` (class 1)
or `2` (class 2). `obdeect stream to-simtel` prints the value for each telescope.

### ATMPROF and synthesised headers

`--atmprof FILE` reads a CORSIKA/sim_telarray atmospheric table (altitude [km], density, thickness,
n-1 per line, `#` comments). It is written as ATMPROF (1216) before the first event, with the
profile number that EVTH names (99 for synthesised headers). A stream with its own ATMPROF block
rejects the option.

A stream without CORSIKA blocks gets synthesised headers. The run end (RUNE) and the input block
(INPUTCFG, one line `NSHOW <n>`) carry the number of events written, so sim_telarray reports the
run as complete:

```text
Run(s) completed as expected after 3 events (3 showers).
```

### Telescope-frame photons (dual-mirror telescopes)

A telescope listed with `mirror_class` 2 may carry telescope-frame (T) photons instead of
camera-frame photons, for example four flashers at the camera firing onto the secondary mirror.
The T frame has its origin at the secondary centre and z along the optical axis towards the sky.
Each photon is an emission point (`x`, `y`, `z`) with direction (`cx`, `cy`, `cz`) and emission
time `t`. `obdeect stream to-simtel` writes one bunch per photon at the plane z = 0 of that frame,
and sim_telarray needs `Bypass_Optics=1`:

```text
photon:  x, y, z = 0.1, 0.2, 0.5 m   c = (0.6, 0, -0.8)   t = 3 ns   observation level 2200 m
bunch:   x, y = 47.5, 20.0 cm   cx, cy = 0.6, 0   time = 5.085 ns   zem = 220050 cm
```

The bunch position is where the photon's line crosses z = 0 (path 0.625 m here, also for a photon
moving away from the plane). The time is `t` plus that path over the speed of light; `zem` is the
observation level plus `z`. Without a `cz` column the photon moves towards the plane. The camera
rotation is not used, and TELPOS z is set to 0 as for camera-frame telescopes.

```console
$ cat cam.txt
1 2 0.0
$ obdeect sources flasher --type secondary --photons 1e6 | obdeect stream to-simtel --camera-map cam.txt | sim_telarray ... -C Bypass_Optics=1 -
obdeect stream to-simtel: telescope 1: telescope frame, mirror class 2, Bypass_Optics=1
```

Errors: T photons for a telescope without map entry or with `mirror_class` 1; a photon with
cz = 0; camera-frame and T photons for the same telescope in one run. The efficiency status of
such a telescope is checked as for camera frame.

### GrISU photon text

`obdeect stream to-grisu` writes ground-frame photons as GrISU/GrOptics photon text. The
atmosphere must be applied upstream and the wavelength assigned, for example by
[`obdeect propagation run`](../propagation/propagation.md):

```console
$ zstd -dc run.corsika.zst | obdeect stream from-corsika \
    | obdeect propagation run --transmission atm.dat --site-altitude 2156 --seed 1 \
        --envelope qe=qe_envelope.dat -i - -o - \
    | obdeect stream to-grisu --telescope-config TLLOC.dat > run.grisu
```

Output (one header block, then `S` and `P` lines):

```text
* HEADF  <-- Start of header flag
* Produced by obdeect stream to-grisu 0.1.0 from a photon stream written by obdeect stream from-corsika 0.1.0
* CORSIKA run number: 5
* Energy range [TeV]: 0.1 10
* Spectral slope: -2
* Primary id: 14
* Primary zenith [deg]: 20
* Primary azimuth [deg, CORSIKA]: 0
* Observation height [m]: 2200
* DATAF  <-- end of header flag
R 1
H 2200
S 1 -50 -100 0 -0.34202014 25000 -1 -1 -1
P -20.5 -11.5 -0.5 +0 +20000 +12.5 400 3 7
```

- `R` is the global QE factor (`--qe`, default 1). It is written, not applied.
- `S E_TeV x y dcos dsin h_first -1 -1 -1`: one line per shower position (per array reuse), with
  core position [m] and primary direction cosines in the GrISU frame. `h_first` is the event key
  `first_interaction_m` [m a.s.l.]. `--event-details` adds `C h_first X_first shower_id`;
  `X_first` is -1 (not carried by the stream) and `shower_id` is the event number.
- `P x y dcx dcy zem t lambda 3 tel`: one line per photon. `x`, `y` [m] are the telescope position
  plus the photon position; `zem` the emission height [m]; `t` the arrival time [ns]; `lambda` the
  wavelength [nm, integer]; `3` the emitter type (unknown); `tel` the GrISU telescope id (1-based).
- A bunch with weight above 1 is an error unless `--expand-bunches` writes `round(w)` copies. A
  weight below 1 is an error: thin the photons upstream.
- Numbers use 9 significant digits. Values below 1e-8 in magnitude are written as 0.

Frames: the CORSIKA ground frame (x north, y west) maps to GrISU by `x_g = -y_c`, `y_g = -x_c`.
Azimuths map by `az_g = (3 pi / 2 - az_c) mod 2 pi`, with `az_c = atan2(cy, cx)` for photons; the
written direction cosines are `sin(ze) cos(az_g)` and `sin(ze) sin(az_g)`, with
`sin(ze) = sqrt(cx^2 + cy^2)`. The primary uses the event zenith and azimuth the same way.

`--telescope-config FILE` reads `* TLLOC id x y z offset px py` lines (x east, y north [m]; other
lines are comments). A stream telescope gets the `id` of the closest line within 0.5 m of its
position (transformed as above). Other telescopes are skipped and reported on stderr. Without the
file, `tel` is the stream telescope id.

In Python: `obdeect.stream.stream_to_grisu(run.ps, run.grisu, telescope_config="TLLOC.dat")`.

### Efficiency status

`obdeect stream to-simtel` refuses a stream whose efficiency status says a factor is still to be
applied when sim_telarray does not apply it in that frame. The error names the factor:

```text
obdeect stream to-simtel: efficiency status of telescope 2 (camera frame): factor
'telescope_transmission_constant' is pending, but sim_telarray does not apply it in this frame;
apply it upstream (state applied); pass --allow-pending telescope_transmission_constant to accept it as it is
```

Use `--allow-pending FACTOR` (repeatable) only when the factor is applied elsewhere in the chain,
for example by the camera module. Factors per frame: [](format.md).

```console
$ obdeect stream to-simtel --camera-map cam.txt --allow-pending telescope_transmission_constant run.ps > sim.simtel
```

In Python: `obdeect.stream.stream_to_simtel(run.ps, sim.simtel, cam, allow_pending=["telescope_transmission_constant"])`.

If the QE envelope was applied upstream (state `envelope`) and the sim_telarray QE table is divided
by that envelope, declare it:

```console
$ obdeect stream to-simtel --camera-map cam.txt --scaled-factor qe run.ps > sim.simtel
```

### Developer checks

`validation/stream/validate_container.sh OUTDIR` is a manual end-to-end check of the codec. In the
`ghcr.io/gammasim/simtools-dev` container it runs CORSIKA, then
`obdeect stream from-corsika | obdeect stream to-simtel`, and compares the CORSIKA EventIO file
with the round-trip file. It then runs sim_telarray on both files with the same configuration and
seed and compares the outputs. The container, CORSIKA and sim_telarray are not part of the
automated tests.

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
consumers and build options: [](../installation.md).

The version comes from the git tag (setuptools-scm, passed by scikit-build-core as
`SKBUILD_PROJECT_VERSION`). Other builds take the newest `v*` tag, or 0.0.0 without one.

## Other languages

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

transform("ground.ps", "camera.ps", trace, stage_id=0x444C0010, frame_out="camera",
          efficiency_status={"mirror": {"state": "applied", "source": "my-model-v1"}})
```

`transform` copies all other records unchanged (records of unknown type are dropped), sets the
frame of the output chunks and telescope records, and appends one `efficiency_status` record per
run. It refuses a factor that is `applied` in any input efficiency status.

The same in C (error handling shortened):

```c
#include <obdeect/stream/obdeect_stream.h>

obdeect_stream_reader *in;  obdeect_stream_writer *out;
obdeect_stream_reader_open("ground.ps", &in);
obdeect_stream_writer_open("camera.ps", "my-tracer", &out);
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

`cpp/tests/data/conformance/` holds `valid.ps` (every record type, plus a record of unknown type
that must be skipped), four invalid streams (bad magic, unknown major version, truncated chunk,
more than 2^16 rows) and `expected.json` with the decoded content and the expected error.
`make_conformance.py` regenerates them with the standard library only. An independent
implementation of [](format.md) should decode `valid.ps` to `expected.json` and reject the others.
