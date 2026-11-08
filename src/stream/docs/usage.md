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
$ ps2simtel --atmprof atmprof1.dat light.ps > sim.simtel  # table for streams without CORSIKA input
$ corsika2ps < run.iact | ps2simtel > sim.simtel       # pipe; "-" reads stdin
$ ps-cat --version                       # ps-cat <package version> (photon-stream format 1.0)
```

The camera map lists the camera-frame telescopes, one per line:

```text
# telescope_id mirror_class(1 single | 2 dual) camera_rotation_deg
2 1 0.0
```

`--atmprof FILE` reads a CORSIKA/sim_telarray atmospheric table (altitude [km], density, thickness, n-1
per line, `#` comments). It is written as ATMPROF (1216) before the first event, with the profile number
that EVTH names (99 for synthesised headers). A stream with its own ATMPROF block rejects the option.

Telescopes not listed are ground-frame (sim_telarray `Bypass_Optics=0`). Listed telescopes with
camera-frame photons are converted to the bypass convention and need `Bypass_Optics=1` (class 1) or
`2` (class 2) in sim_telarray; `ps2simtel` prints the value for each. Listed dual-mirror telescopes
may instead carry telescope-frame photons (next section).

### Telescope-frame photons (dual-mirror telescopes)

A dual-mirror telescope listed with `mirror_class` 2 may carry telescope-frame (T) photons
instead of camera-frame photons, for example four flashers at the camera firing onto the
secondary mirror. The T frame has its origin at the secondary centre and z along the optical
axis towards the sky; each photon is an emission point (`x`, `y`, `z`) with direction
(`cx`, `cy`, `cz`) and emission time `t`. `ps2simtel` writes one bunch per photon at the
plane z = 0 of that frame, and sim_telarray needs `Bypass_Optics=1`:

```text
photon:  x, y, z = 0.1, 0.2, 0.5 m   c = (0.6, 0, -0.8)   t = 3 ns   observation level 2200 m
bunch:   x, y = 47.5, 20.0 cm   cx, cy = 0.6, 0   time = 5.085 ns   zem = 220050 cm
```

The bunch position is where the photon's line crosses z = 0 (path 0.625 m here, also for a
photon moving away from the plane), the time is `t` plus that path length over the speed of
light, and `zem` is the observation level plus `z`. Without a `cz` column the photon moves
towards the plane. The camera rotation of the map is not used, and TELPOS z is set to 0 as for
camera-frame telescopes. `ps2simtel` reports the value when it sees the first photons:

```console
$ cat cam.txt
1 2 0.0
$ dlamb-flasher --type secondary --photons 1e6 | ps2simtel --camera-map cam.txt | sim_telarray ... -C Bypass_Optics=1 -
ps2simtel: telescope 1: telescope frame, mirror class 2, Bypass_Optics=1
```

Errors: T photons for a telescope without map entry or with `mirror_class` 1, a photon with
cz = 0, and camera-frame and T photons for the same telescope in one run.
The efficiency ledger of such a telescope is checked as for camera frame.

### Efficiency ledger

`ps2simtel` refuses a stream whose ledger says a factor is still to be applied and
sim_telarray does not apply it in that frame. The error names the factor:

```text
ps2simtel: efficiency ledger of telescope 2 (camera frame): factor
'telescope_transmission_constant' is pending, but sim_telarray does not apply it in this frame;
apply it upstream (state applied); pass --allow-pending telescope_transmission_constant to accept it as it is
```

Use `--allow-pending FACTOR` (repeatable) only when the factor is applied elsewhere
in the chain, for example by the camera module. Which factors sim_telarray applies in
each frame is listed in docs/format.md ("Efficiency ledger").

```console
$ ps2simtel --camera-map cam.txt --allow-pending telescope_transmission_constant run.ps > sim.simtel
```

In Python: `gstrickd.stream_to_simtel(run.ps, sim.simtel, cam, allow_pending=["telescope_transmission_constant"])`.

When the QE envelope was applied upstream (ledger state `envelope`) and the
sim_telarray QE table is divided by that envelope, declare it:

```console
$ ps2simtel --camera-map cam.txt --scaled-factor qe run.ps > sim.simtel
```

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
