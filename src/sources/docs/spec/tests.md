# Tests and required API

The tests are the acceptance criteria. C++ tests are `cpp/tests/test_<name>.cpp`, one
executable and CTest each, built with `check.hpp` (`CHECK`, `CHECK_NEAR`, `run_checks`);
helpers `program_stream.hpp` (`read_stream` -> `StreamContents{types, metadata, rows}`,
`Row`, `count_type`) and `sky_stream_check.hpp` (`read_stream` -> `StreamSummary{types,
photons per telescope, chunks, events, run metadata, ended}`, `write_to_string`) read streams
back with `gstrickd::Reader`. All are kept unchanged. Python tests: `python -m pytest python/tests`.

## Required C++ API (as used by the tests)

| Header | Symbols (namespace `dlamb` unless stated) |
| --- | --- |
| `density.hpp` | `Uniforms`, `pi`, `two_pi`, `normal_quantile`, `LinearDensity(x, y)` (`integral`, `maximum`, `at`, `cdf`, `sample`), `Fixed`, `UniformRange`, `GaussianShape{mean_value, sigma}`, `Exponential{mean_value, tau}`, `ComplexPulse{flat_width, gauss_rms, decay, delay}`, `InverseSquare{low, high}` (`mean()`, `variance()`), `Tabulated{density}`, `Isotropic`, `FixedAngle{u0}`, `FilledCone{u0}`, `Lambertian`, `TabulatedAngular{density}` (each with `cdf`, `quantile`), variants `Spectrum`, `Pulse`, `Angular`, `angle_quantile`, `gaussian_angular_table`, `rayleigh_angular_table` |
| `source.hpp` | `Vec3` helpers (`vec3.hpp`: `dot`, `norm`, `normalised`), `speed_of_light_m_per_ns`, `Emitter`, `Star`, `DiscTarget`, `SphereTarget`, `Source`, `Photon`, `source_columns`, `max_photon_count`, `validate`, `expected_count`, `photon_count(source, seed, event, telescope)`, `emit(source, seed, event, telescope, first_index, rows, PhotonChunk&)`, `sample_photon` |
| `streams.hpp` | `stage_sources`, `photon_stream_format_major`, `Stream` |
| `track.hpp` | `Track`, `TrackArc`, `track_arc(track, s)`, `track_count`, `emit_track`, `expected_track_weight`, `validate(Track)` |
| `emission_stream.hpp` | `source_seed`, `StreamSummary{photons, weight, expected}`, `RunWriter`, `RunInfo`, `EventInfo`, `TelescopeSphere` |
| `background.hpp` | `BackgroundSource`, `sample_background`, `emit_background`, `background_count` |
| `sky_sources.hpp` | `SkyRun`, `SkyTelescope`, `SourceEmission`, `SkyPlan`, `write_plan`, `shared_telescope_key`, `direction_from_altitude`, `rotate_y`, `rotate_z`, `fiducial_super_radius` |
| `cli.hpp` | `dlamb::cli::parse_number`, `parse_integer` |
| `programs/stream_output.hpp` | `dlamb::programs`: `stage_programs`, `Emission`, `Event`, `RunSetup`, `derive_seed`, `write_stream(out, RunSetup, events)`, `parse_*` helpers |
| `programs/flasher.hpp` | `programs::flasher`: `Type`, `Options` (fields in the flasher page), `parse_type`, `require_type`, `observation_level_m`, `build_sources`, `write_stream(out, Options)` |
| `programs/led.hpp` | `programs::led`: `Options`, `Camera`, `read_camera(istream, keep_slaves)`, `build_sources`, `write_stream`, `Event::labels` |
| `programs/point-source.hpp` | `programs::XyzlsOptions` (`events`, `photons`, `bunch_size`, `z_cm`, `altitude_m`, `angular`, `pulse`, `spectrum`, ...), `TelposCm`, `plan_point_source`, `write_point_source` |
| `programs/beam-scan.hpp` | `programs::XylsOptions` (`radius_cm`, `x_cm`, `y_cm`, `height_cm`, `photons_per_position`, `bunch_size`), `plan_beam_scan`, `write_beam_scan` |
| `programs/laser-beam.hpp` | `programs::LsBeamOptions`, `BeamGeometry{start_m, direction, path_m}`, `laser_beam_geometry`, `plan_laser_beam`, `write_laser_beam` |
| `programs/night-sky.hpp` | `programs::NsbOptions`, `plan_night_sky`, `write_night_sky`; plan emissions expose `count(event)`, `fill(event, first, rows, chunk)` |
| `programs/muon-track.hpp` | `muon_track::Options`, `IndexProfile`, `tracks`, `validate`, `alpha_fine`, `event_count`, `primary_id_from_name`, `mass_charge`, `read_index_profile`, `parse_command_line`, `write` |
| `tools/sky_args.hpp` | `dlamb::cli`: `to_count`, `to_double`, `parse_spectrum`, `parse_pulse(text, option, bare_is_gauss)`, `parse_angular`, `is_number` |
| `tools/telpos_file.hpp` | `dlamb::cli::read_telpos(istream, option, corsika)` |

Struct field names and defaults are those of the corresponding tool page (`Options` of a tool:
the option table with units as stated; `XyzlsOptions`: `run_number`, `events`,
`telcam_radius_cm`, `x_cm`, `y_cm`, `z_cm`, `photons` (vector), `bunch_size`, `altitude_m`,
`atmosphere`, `direction`, `spectrum`, `pulse`, `angular`, `seed`, `telescopes`;
`LsBeamOptions`: `events`, `step_m`, `bunch_size`, `primary_photons_per_m`, `spectrum`,
`pulse`, `angular` (string), `run_number`, `altitude_m`, `atmosphere`, `laser_position_m`,
`laser_theta_deg`, `laser_phi_deg`, `telescope_theta_deg`, `telescope_phi_deg`,
`telescope_radius_m`, `fov_deg`, `seed`; `NsbOptions`: `telescope_radius_m`, `viewcone_deg`,
`altitude_m`, `atmosphere`, `photons`, `bunch_size`, `events`, `run_number`, `spectrum`,
`pulse_width_ns`, `pulse_delay_ns`, `seed`; `flasher::Options`: `type`, `frame`, `photons`
(vector), `bunch_size`, `spectrum`, `lightpulse`, `angular` (strings), `run_number`,
`altitude_m` (optional), `atmosphere`, `seed`, `events` (vector<int>), `x_cm`, `y_cm`,
`distance_cm`, `camera_radius_cm`, `flasher_xy_cm`, `flasher_depth_cm`, `inclination_deg`,
`mirror_camera_distance_cm`, `delays` (array of 5), `fire`; `led::Options`: `camera`, `slaves`,
`use`, `trigger_random`, `trigger_all`, `require`, `random`, `events`, `photons`, `bunch_size`,
`camera_radius_cm`, `spectrum`, `lightpulse`, `angular` ("parallel"), `run_number`,
`altitude_m`, `atmosphere`, `seed`; muon `Options`: see tool page, fields `run_number`, `seed`,
`events`, `observation_level_m`, `telescope_radius_m`, `impact_x_m`, `impact_y_m`, `rscat_m`,
`nscat`, `h_start_m`, `theta_deg`, `phi_deg`, `index_minus_one`, `index_profile`, `bunch_size`,
`wavelength_min_nm`, `wavelength_max_nm`, `charge`, `beta`, `momentum_gev` (optional),
`mass_gev`, `primary_id`, `multiple_scattering`).

## C++ tests

| Test | What it checks |
| --- | --- |
| `test_version` | `version_major/minor/patch` constants exist |
| `test_density` | `normal_quantile` references; `LinearDensity` integral, cdf, sample, rejection of unsorted nodes; moments of every 1-D shape and angular profile at four sigma; Gaussian table values; widths 1e-200, 0, -1, inf rejected |
| `test_source` | star plane wave and point source geometry and times; disc target weights, candidate count, angular cut; chunk-size, order, seed, event independence; chunk start at 2^20; rows above 65536 rejected; unit directions; counts (floor plus Bernoulli, cone target count); sphere target hits; every `validate` rule and bound |
| `test_gstrickd_contract` | stage ids, frame codes, column names, chunk limit, format major |
| `test_sky_args` | `to_count`, fixed-wavelength range, Cherenkov syntax, pulse names and arity, angular names, exponent numbers |
| `test_flasher` | dish counts, moments, zem, bunches, event numbering, camera frame; secondary positions, counts, delays, fire patterns; type rules; non-finite inputs |
| `test_led` | see the led page |
| `test_point-source` | see the point-source page; telescope file reader |
| `test_beam-scan` | position count by repeated addition, layout, bunches, rejections |
| `test_laser-beam` | see the laser-beam page |
| `test_night-sky` | see the night-sky page |
| `test_muon-track` | see the muon-track page |

## Python tests

| Test file | What it checks |
| --- | --- |
| `test_source.py` | float32 columns, unit directions, disc weights, plane-wave star, finite-distance star, angular mean cosines (cone, Lambertian, isotropic), Gaussian and exponential pulse moments, tables, Gaussian angular width, chunk independence, `iter_chunks`, counts, sphere target, every rejection listed in the python page, cone angle range, table angle range, representable Gauss width, counts above 2^53, complex pulse widths |
| `test_gstrickd_contract.py` | Python constants equal the compiled gstrickd values (and the installed package when importable) |
| `test_plot.py` | three panels with the stated labels (skipped without matplotlib) |
| `test_import.py` | package and extension import; `core_version()` has three entries |

## Tools

No test runs the executables; the tool `main` functions only parse options and call the library
functions above. Tool behaviour that no test covers (option parsing, messages, `--list-only`)
is specified on the tool pages and checked by hand.

## `dlamb/version.hpp`

`dlamb::version_major`, `version_minor`, `version_patch` (non-negative integers, matching the
project version); the Python `core_version()` returns them as a tuple.

## Not specified at code level

- Exact Philox outputs are defined by gstrickd (`uniform4`), not here; the stream
  assignments and keys are specified in [conventions](conventions.md).
- The numerically stable quadratic solver of the laser-beam and the normal quantile
  approximation are specified by their accuracy only.
- The night-sky default table (`Benn_LaPalma_sky_converted.lis`) is external data.
- Partial output after a late source-validation error, and usage text wording, are not
  contractual.
