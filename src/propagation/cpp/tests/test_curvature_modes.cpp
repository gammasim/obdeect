// Curvature modes of the atmosphere stage (propagator): input checks and rejections, efficiency
// status keys, the CORSIKA CURVED extinction, the spherical sky geometry (limit R -> infinity and
// the Earth-radius correction), the spherical refraction and air delay with a profile, the
// dropped photons above the telescope level, and independence of chunk size and row order in the
// curved modes.
#include "check.hpp"
#include "stream_fixtures.hpp"

#include "obdeect/propagation/curved_refraction.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/propagator.hpp"
#include "obdeect/propagation/refraction.hpp"
#include "obdeect/propagation/scattering.hpp"
#include "obdeect/propagation/sky.hpp"
#include "obdeect/propagation/timing.hpp"
#include "obdeect/propagation/transmission.hpp"

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <variant>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using namespace obdeect_test;

using obdeect::propagation::Curvature;
using obdeect::propagation::OutputMode;
using obdeect::propagation::PropagatorConfig;
using obdeect::stream::Column;
using obdeect::stream::Frame;
using obdeect::stream::RecordType;

constexpr double site_m = 2000.0;
constexpr double kEarthRadius = obdeect::stream::earth_radius_m;
constexpr std::int64_t no_key = -1; // run key curvature_applied absent
const double kPi = obdeect::stream::pi;

obdeect::propagation::TransmissionTable real_table() {
  return obdeect::propagation::TransmissionTable::from_file(
      data_file("atm_trans_2000_1_10_0_0_2000.dat"), site_m,
      obdeect::propagation::WavelengthInterpolation::floor);
}

// Two telescopes: at the array reference (z = 0) and 3 km east, 2.5 km south of it.
// curvature_applied: 0 or 1 as written, no_key: the key is absent.
obdeect::stream::Metadata run_keys(std::int64_t curvature_applied) {
  obdeect::stream::Metadata m;
  m["run_number"] = std::int64_t{1};
  m["observation_level_m"] = site_m;
  m["telescope_x_m"] = std::vector<double>{100.0, 3000.0};
  m["telescope_y_m"] = std::vector<double>{0.0, -2500.0};
  m["telescope_z_m"] = std::vector<double>{0.0, 0.0};
  if (curvature_applied != no_key) {
    m["curvature_applied"] = curvature_applied;
  }
  m["wavelength_min_nm"] = 250.0;
  m["wavelength_max_nm"] = 700.0;
  return m;
}

// Frame G photon. The wavelength is given, so no wavelength is drawn.
struct GPhoton {
  double x{}, y{}, cx{}, cy{}, zem{}, wavelength{}, w{1.0};
};

// Sky-frame photon (cz = -sqrt(1 - cx^2 - cy^2) unless upward).
struct SkyPhoton {
  double x{}, y{}, z{}, cx{}, cy{}, w{1.0};
  bool upward{};
  double wavelength{}; // [nm] as input column when > 0 (else drawn from the photon index)
};

// Downward photons at zenith angles from 0 to 85 degrees. x is the photon id (kept in the output).
std::vector<GPhoton> ground_photons(std::size_t n) {
  std::vector<GPhoton> photons;
  for (std::size_t i = 0; i < n; ++i) {
    const double zenith =
        85.0 * kPi / 180.0 * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
    const double phi = 0.7 * static_cast<double>(i);
    GPhoton photon;
    photon.x = static_cast<double>(i);
    photon.y = 100.0 * std::sin(0.37 * static_cast<double>(i));
    photon.cx = std::sin(zenith) * std::cos(phi);
    photon.cy = std::sin(zenith) * std::sin(phi);
    photon.zem = 2000.0 + 10000.0 * static_cast<double>(i % 7) / 6.0;
    photon.wavelength = 300.0 + static_cast<double>(i % 400);
    photons.push_back(photon);
  }
  return photons;
}

// Sky-frame photons 2 to 17 km above the observation level; x is the lateral position.
std::vector<SkyPhoton> sky_photons(std::size_t n) {
  std::vector<SkyPhoton> photons;
  for (std::size_t i = 0; i < n; ++i) {
    const auto k = static_cast<double>(i);
    SkyPhoton photon;
    photon.x = 1500.0 * std::sin(0.1 * k);
    photon.y = 1500.0 * std::cos(0.07 * k);
    photon.z = 3000.0 + 4000.0 * static_cast<double>(i % 5);
    photon.cx = 0.3 * std::sin(0.01 * k);
    photon.cy = 0.3 * std::cos(0.013 * k);
    photon.w = 1.0 + static_cast<double>(i % 3);
    photon.upward = i % 17 == 0;
    photons.push_back(photon);
  }
  return photons;
}

obdeect::stream::PhotonChunk ground_chunk(std::uint32_t telescope, std::size_t first,
                                          const std::vector<GPhoton> &photons, std::size_t rows) {
  obdeect::stream::PhotonChunk chunk;
  chunk.telescope_id = telescope;
  chunk.frame = Frame::ground;
  chunk.first_index = first;
  const std::uint32_t mask =
      obdeect::stream::column_bit(Column::x) | obdeect::stream::column_bit(Column::y) |
      obdeect::stream::column_bit(Column::cx) | obdeect::stream::column_bit(Column::cy) |
      obdeect::stream::column_bit(Column::t) | obdeect::stream::column_bit(Column::zem) |
      obdeect::stream::column_bit(Column::w) | obdeect::stream::column_bit(Column::wavelength);
  chunk.resize(rows, mask);
  for (std::size_t r = 0; r < rows; ++r) {
    const GPhoton &photon = photons[first + r];
    chunk[Column::x][r] = static_cast<float>(photon.x);
    chunk[Column::y][r] = static_cast<float>(photon.y);
    chunk[Column::cx][r] = static_cast<float>(photon.cx);
    chunk[Column::cy][r] = static_cast<float>(photon.cy);
    chunk[Column::t][r] = 0.0F;
    chunk[Column::zem][r] = static_cast<float>(photon.zem);
    chunk[Column::w][r] = static_cast<float>(photon.w);
    chunk[Column::wavelength][r] = static_cast<float>(photon.wavelength);
  }
  return chunk;
}

obdeect::stream::PhotonChunk sky_chunk(std::uint32_t telescope, std::size_t first,
                                       const std::vector<SkyPhoton> &photons, std::size_t rows) {
  obdeect::stream::PhotonChunk chunk;
  chunk.telescope_id = telescope;
  chunk.frame = Frame::sky;
  chunk.first_index = first;
  // A wavelength column is written when the first photon of the chunk has one (all do in a stream).
  const bool has_wavelength = photons[first].wavelength > 0.0;
  std::uint32_t mask =
      obdeect::stream::column_bit(Column::x) | obdeect::stream::column_bit(Column::y) |
      obdeect::stream::column_bit(Column::z) | obdeect::stream::column_bit(Column::cx) |
      obdeect::stream::column_bit(Column::cy) | obdeect::stream::column_bit(Column::cz) |
      obdeect::stream::column_bit(Column::t) | obdeect::stream::column_bit(Column::w);
  if (has_wavelength) {
    mask |= obdeect::stream::column_bit(Column::wavelength);
  }
  chunk.resize(rows, mask);
  for (std::size_t r = 0; r < rows; ++r) {
    const SkyPhoton &photon = photons[first + r];
    const double magnitude = std::sqrt(1.0 - photon.cx * photon.cx - photon.cy * photon.cy);
    chunk[Column::x][r] = static_cast<float>(photon.x);
    chunk[Column::y][r] = static_cast<float>(photon.y);
    chunk[Column::z][r] = static_cast<float>(photon.z);
    chunk[Column::cx][r] = static_cast<float>(photon.cx);
    chunk[Column::cy][r] = static_cast<float>(photon.cy);
    chunk[Column::cz][r] = static_cast<float>(photon.upward ? magnitude : -magnitude);
    chunk[Column::t][r] = 0.0F;
    chunk[Column::w][r] = static_cast<float>(photon.w);
    if (has_wavelength) {
      chunk[Column::wavelength][r] = static_cast<float>(photon.wavelength);
    }
  }
  return chunk;
}

// Stream of one event: every photon goes to telescopes 1 and 2, in chunks of chunk_rows.
std::string ground_stream(const obdeect::stream::Metadata &run, const std::vector<GPhoton> &photons,
                          std::size_t chunk_rows) {
  std::ostringstream out;
  {
    obdeect::stream::Writer writer(out, "test");
    writer.write_metadata(RecordType::run, run);
    writer.write_metadata(RecordType::event, event_keys());
    for (const std::uint32_t telescope : {1U, 2U}) {
      writer.write_metadata(RecordType::telescope, telescope_keys(telescope, Frame::ground));
      for (std::size_t first = 0; first < photons.size(); first += chunk_rows) {
        const std::size_t rows = std::min(chunk_rows, photons.size() - first);
        writer.write_chunk(ground_chunk(telescope, first, photons, rows));
      }
      writer.write_metadata(RecordType::telescope_end, telescope_end_keys(telescope));
    }
    writer.write_metadata(RecordType::event_end, {});
    writer.write_metadata(RecordType::run_end, {});
    writer.finish();
  }
  return out.str();
}

std::string sky_stream(const obdeect::stream::Metadata &run, const std::vector<SkyPhoton> &photons,
                       std::size_t chunk_rows) {
  std::ostringstream out;
  {
    obdeect::stream::Writer writer(out, "test");
    writer.write_metadata(RecordType::run, run);
    writer.write_metadata(RecordType::event, event_keys());
    for (const std::uint32_t telescope : {1U, 2U}) {
      writer.write_metadata(RecordType::telescope, telescope_keys(telescope, Frame::sky));
      for (std::size_t first = 0; first < photons.size(); first += chunk_rows) {
        const std::size_t rows = std::min(chunk_rows, photons.size() - first);
        writer.write_chunk(sky_chunk(telescope, first, photons, rows));
      }
      writer.write_metadata(RecordType::telescope_end, telescope_end_keys(telescope));
    }
    writer.write_metadata(RecordType::event_end, {});
    writer.write_metadata(RecordType::run_end, {});
    writer.finish();
  }
  return out.str();
}

// Column c keyed by the column x of the same rows (telescope of the output).
std::map<float, float> value_by_x(const std::vector<obdeect::stream::Record> &records,
                                  std::uint32_t telescope, Column c) {
  std::map<float, float> values;
  for (const auto &record : records) {
    if (record.type != RecordType::photons || record.chunk.telescope_id != telescope ||
        !record.chunk.has(c)) {
      continue;
    }
    const auto &x = record.chunk[Column::x];
    const auto &value = record.chunk[c];
    for (std::size_t row = 0; row < value.size(); ++row) {
      values[x[row]] = value[row];
    }
  }
  return values;
}

// x of the photons of one telescope in the order of the output (sorted: thinned selections).
std::vector<float> sorted_x(const std::vector<obdeect::stream::Record> &records,
                            std::uint32_t telescope) {
  auto x = photon_column(records, telescope, Column::x);
  std::sort(x.begin(), x.end());
  return x;
}

const obdeect::stream::Metadata *
first_record_of(const std::vector<obdeect::stream::Record> &records, RecordType type) {
  for (const auto &record : records) {
    if (record.type == type) {
      return &record.metadata;
    }
  }
  return nullptr;
}

bool has_key(const obdeect::stream::Metadata &m, const std::string &key) {
  return m.find(key) != m.end();
}

std::string text_of(const obdeect::stream::Metadata &m, const std::string &key) {
  const auto it = m.find(key);
  if (it == m.end() || !std::holds_alternative<std::string>(it->second)) {
    return "<missing>";
  }
  return std::get<std::string>(it->second);
}

double number_of(const obdeect::stream::Metadata &m, const std::string &key) {
  const auto it = m.find(key);
  if (it == m.end()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (const auto *value = std::get_if<double>(&it->second)) {
    return *value;
  }
  return static_cast<double>(std::get<std::int64_t>(it->second));
}

struct Result {
  std::string bytes;
  obdeect::propagation::PropagationSummary summary;
};

Result propagate(const std::string &input, const PropagatorConfig &config) {
  std::istringstream in(input);
  std::ostringstream out;
  Result result;
  result.summary = obdeect::propagation::propagate_stream(in, out, config);
  result.bytes = out.str();
  return result;
}

// Weighted output with the acceptance column, the default for these tests.
PropagatorConfig config_of(Curvature curvature, double earth_radius_m) {
  PropagatorConfig config;
  config.transmission.emplace(real_table());
  config.transmission_source = "test table";
  config.site_altitude_m = site_m;
  config.seed = 20261114;
  config.mode = OutputMode::weighted;
  config.store_p = true;
  config.curvature = curvature;
  config.earth_radius_m = earth_radius_m;
  return config;
}

// The atmospheric profile of test_propagator.cpp (same site altitude).
obdeect::propagation::AtmosphereProfile ecmwf_profile() {
  return obdeect::propagation::AtmosphereProfile::from_file(
      data_file("atmprof_ecmwf_north_winter_fixed.dat"));
}

// Threads of the propagator and of its curved tables: the output does not depend on the count.
unsigned build_threads() { return std::max(2U, std::thread::hardware_concurrency()); }

// Weighted output with the acceptance column and the profile (spherical or CORSIKA curved runs).
PropagatorConfig profile_config(Curvature curvature, double earth_radius_m) {
  PropagatorConfig config = config_of(curvature, earth_radius_m);
  config.profile.emplace(ecmwf_profile());
  config.profile_source = "ecmwf";
  config.threads = build_threads();
  return config;
}

// Run keys with telescopes at (x_m, y_m, 0) and the key refraction_applied (profile runs).
obdeect::stream::Metadata lateral_run_keys(std::int64_t curvature_applied,
                                           std::int64_t refraction_applied,
                                           const std::vector<double> &x_m,
                                           const std::vector<double> &y_m) {
  obdeect::stream::Metadata m = run_keys(curvature_applied);
  m["telescope_x_m"] = x_m;
  m["telescope_y_m"] = y_m;
  m["telescope_z_m"] = std::vector<double>(x_m.size(), 0.0);
  m["refraction_applied"] = refraction_applied;
  return m;
}

// As above with the telescope altitudes z_m (m above the observation level).
obdeect::stream::Metadata lateral_run_keys(std::int64_t curvature_applied,
                                           std::int64_t refraction_applied,
                                           const std::vector<double> &x_m,
                                           const std::vector<double> &y_m,
                                           const std::vector<double> &z_m) {
  obdeect::stream::Metadata m = lateral_run_keys(curvature_applied, refraction_applied, x_m, y_m);
  m["telescope_z_m"] = z_m;
  return m;
}

// One block of a sky stream: a run record with its own keys, one event with the photons sent to
// telescopes 1 to `telescopes`, in chunks of chunk_rows rows.
struct SkyBlock {
  obdeect::stream::Metadata run;
  std::uint32_t event;
  std::vector<SkyPhoton> photons;
  std::size_t chunk_rows;
  std::uint32_t telescopes;
};

std::string sky_blocks(const std::vector<SkyBlock> &blocks) {
  std::ostringstream out;
  {
    obdeect::stream::Writer writer(out, "test");
    for (const auto &block : blocks) {
      writer.write_metadata(RecordType::run, block.run);
      obdeect::stream::Metadata event;
      event["event_number"] = static_cast<std::int64_t>(block.event);
      writer.write_metadata(RecordType::event, event);
      for (std::uint32_t telescope = 1; telescope <= block.telescopes; ++telescope) {
        writer.write_metadata(RecordType::telescope, telescope_keys(telescope, Frame::sky));
        for (std::size_t first = 0; first < block.photons.size(); first += block.chunk_rows) {
          const std::size_t rows = std::min(block.chunk_rows, block.photons.size() - first);
          writer.write_chunk(sky_chunk(telescope, first, block.photons, rows));
        }
        writer.write_metadata(RecordType::telescope_end, telescope_end_keys(telescope));
      }
      writer.write_metadata(RecordType::event_end, {});
      writer.write_metadata(RecordType::run_end, {});
    }
    writer.finish();
  }
  return out.str();
}

// The records of the output, split before each event record: [0] holds the records before the
// first event, [n] the records of event n (event numbers 1, 2, ... in the order of the input).
std::vector<std::vector<obdeect::stream::Record>>
split_events(const std::vector<obdeect::stream::Record> &records) {
  std::vector<std::vector<obdeect::stream::Record>> events(1);
  for (const auto &record : records) {
    if (record.type == RecordType::event) {
      events.emplace_back();
    }
    events.back().push_back(record);
  }
  return events;
}

// Sky-frame photon at (x, y, z) above the observation level, zenith angle zenith_deg of its
// direction of motion (downwards) and azimuth phi_rad.
SkyPhoton sky_at(double x, double y, double z_m, double zenith_deg, double phi_rad = 0.0) {
  const double zenith = zenith_deg * kPi / 180.0;
  SkyPhoton photon;
  photon.x = x;
  photon.y = y;
  photon.z = z_m;
  photon.cx = std::sin(zenith) * std::cos(phi_rad);
  photon.cy = std::sin(zenith) * std::sin(phi_rad);
  return photon;
}

// Input checks and config checks of the curvature modes.
void test_rejections() {
  const auto ground_in = [](std::int64_t key) {
    return ground_stream(run_keys(key), ground_photons(10), 10);
  };
  const auto sky_in = [](std::int64_t key) {
    return sky_stream(run_keys(key), sky_photons(10), 10);
  };
  const double R = kEarthRadius;

  // none: CORSIKA CURVED input is refused; a missing key and 0 are accepted.
  CHECK(throws_with([&] { propagate(ground_in(1), config_of(Curvature::none, R)); },
                    "CORSIKA CURVED"));
  CHECK(accepts([&] { propagate(ground_in(no_key), config_of(Curvature::none, R)); }));
  CHECK(accepts([&] { propagate(ground_in(0), config_of(Curvature::none, R)); }));

  // corsika: the key must be present and 1.
  CHECK(throws_with([&] { propagate(ground_in(no_key), config_of(Curvature::corsika, R)); },
                    "curvature_applied"));
  CHECK(throws_with([&] { propagate(ground_in(0), config_of(Curvature::corsika, R)); }, "CORSIKA"));
  CHECK(accepts([&] { propagate(ground_in(1), config_of(Curvature::corsika, R)); }));
  // corsika: sky-frame telescopes have no CORSIKA CURVED geometry.
  CHECK(throws_with([&] { propagate(sky_in(1), config_of(Curvature::corsika, R)); }, "sky"));

  // spherical: the key must be present and 0.
  CHECK(throws_with([&] { propagate(ground_in(no_key), config_of(Curvature::spherical, R)); },
                    "curvature_applied"));
  CHECK(throws_with([&] { propagate(ground_in(1), config_of(Curvature::spherical, R)); },
                    "already applied"));
  CHECK(accepts([&] { propagate(sky_in(0), config_of(Curvature::spherical, R)); }));

  // Scattering and a profile are not combined with the curved modes (checked at construction).
  for (const Curvature curvature : {Curvature::corsika, Curvature::spherical}) {
    auto scattering = config_of(curvature, R);
    scattering.scattering.enabled = true;
    scattering.scattering.points_per_photon = 8;
    scattering.scattering_table.emplace(
        obdeect::propagation::ScatteringTable::from_file(data_file("scattering_exponential.dat")));
    CHECK(throws_with(
        [&] {
          obdeect::propagation::Propagator propagator(scattering);
          (void)propagator;
        },
        "curved atmosphere cannot be combined"));
  }
  // A profile is accepted in the spherical and the CORSIKA curved modes (refraction and air delay
  // of the spherical mode: test_profile_* below). A spherical run without refraction_applied is
  // refused when the run record is read.
  for (const Curvature curvature : {Curvature::corsika, Curvature::spherical}) {
    auto with_profile = profile_config(curvature, R);
    CHECK(accepts([&] {
      obdeect::propagation::Propagator propagator(with_profile);
      (void)propagator;
    }));
  }
  CHECK(throws_with(
      [&] {
        propagate(sky_stream(run_keys(0), sky_photons(10), 10),
                  profile_config(Curvature::spherical, R));
      },
      "profile"));

  // Earth radius: finite and above 1e6 m, in every mode.
  for (const Curvature curvature : {Curvature::none, Curvature::corsika, Curvature::spherical}) {
    for (const double radius : {0.0, -1.0, 1e6, std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::infinity()}) {
      CHECK(throws_with(
          [&] {
            obdeect::propagation::Propagator propagator(config_of(curvature, radius));
            (void)propagator;
          },
          "earth radius must be finite"));
    }
  }

  // Mode names.
  CHECK(obdeect::propagation::parse_curvature("none") == Curvature::none);
  CHECK(obdeect::propagation::parse_curvature("corsika") == Curvature::corsika);
  CHECK(obdeect::propagation::parse_curvature("spherical") == Curvature::spherical);
  for (const char *name : {"", "Corsika", "curved", "flat"}) {
    CHECK(throws_with([&] { obdeect::propagation::parse_curvature(name); }, "curvature"));
  }

  // Geometry only: the curved geometry is allowed (no extinction).
  auto geometry_spherical = config_of(Curvature::spherical, R);
  geometry_spherical.geometry_only = true;
  geometry_spherical.transmission.reset();
  geometry_spherical.store_p = false;
  CHECK(accepts([&] { propagate(sky_in(0), geometry_spherical); }));
  auto geometry_corsika = config_of(Curvature::corsika, R);
  geometry_corsika.geometry_only = true;
  geometry_corsika.transmission.reset();
  geometry_corsika.store_p = false;
  CHECK(accepts([&] { propagate(ground_in(1), geometry_corsika); }));
}

// Efficiency status keys per mode; the run record of the spherical output says the curvature is
// applied.
void test_status_and_run_key() {
  const double R = 6.0e6;
  const auto corsika = records_of(propagate(ground_stream(run_keys(1), ground_photons(20), 20),
                                            config_of(Curvature::corsika, R))
                                      .bytes);
  const auto *status = first_record_of(corsika, RecordType::efficiency_status);
  CHECK(status != nullptr);
  if (status != nullptr) {
    CHECK(text_of(*status, "atmosphere.geometry") == "corsika-curved");
    CHECK(number_of(*status, "atmosphere.earth_radius_m") == R);
  }

  const auto spherical = records_of(
      propagate(sky_stream(run_keys(0), sky_photons(20), 20), config_of(Curvature::spherical, R))
          .bytes);
  status = first_record_of(spherical, RecordType::efficiency_status);
  CHECK(status != nullptr);
  if (status != nullptr) {
    CHECK(text_of(*status, "atmosphere.geometry") == "spherical");
    CHECK(number_of(*status, "atmosphere.earth_radius_m") == R);
  }
  const auto *run = first_record_of(spherical, RecordType::run);
  CHECK(run != nullptr);
  if (run != nullptr) {
    CHECK(std::get<std::int64_t>(run->at("curvature_applied")) == 1);
  }

  // none: no curvature keys, and no curvature_applied key in the run record.
  const auto none = records_of(propagate(ground_stream(run_keys(no_key), ground_photons(20), 20),
                                         config_of(Curvature::none, R))
                                   .bytes);
  status = first_record_of(none, RecordType::efficiency_status);
  CHECK(status != nullptr);
  if (status != nullptr) {
    CHECK(!has_key(*status, "atmosphere.geometry"));
    CHECK(!has_key(*status, "atmosphere.earth_radius_m"));
  }
  run = first_record_of(none, RecordType::run);
  CHECK(run != nullptr);
  if (run != nullptr) {
    CHECK(!has_key(*run, "curvature_applied"));
  }
}

// One photon at zenith 80 degrees (and one at zenith 0): the curved acceptance is the table's
// transmission_curved, larger than the plane-parallel value at 80 degrees, equal at zenith 0.
void test_corsika_extinction() {
  const double R = kEarthRadius;
  const auto table = real_table();
  const double zenith = 80.0 * kPi / 180.0;
  GPhoton slanted;
  slanted.cx = std::sin(zenith);
  slanted.zem = 10000.0;
  slanted.wavelength = 400.0;
  GPhoton vertical = slanted;
  vertical.cx = 0.0;

  const auto acceptance = [&](std::int64_t key, Curvature curvature, const GPhoton &photon) {
    const auto records = records_of(
        propagate(ground_stream(run_keys(key), {photon}, 1), config_of(curvature, R)).bytes);
    const auto p = photon_column(records, 1, Column::p);
    CHECK(p.size() == 1);
    return p.empty() ? -1.0 : static_cast<double>(p[0]);
  };

  const double p_curved = acceptance(1, Curvature::corsika, slanted);
  const double p_flat = acceptance(no_key, Curvature::none, slanted);
  const double expected =
      table.transmission_curved<double>(400.0, 10000.0, std::cos(zenith), site_m, 0.0, R);
  CHECK_NEAR(p_curved, expected, 1e-5);
  CHECK(p_curved > p_flat);
  std::fprintf(stderr, "corsika zenith 80 deg: curved %.6f, plane-parallel %.6f\n", p_curved,
               p_flat);

  CHECK_NEAR(acceptance(1, Curvature::corsika, vertical),
             acceptance(no_key, Curvature::none, vertical), 1e-5);
}

// Spherical geometry with R = 1e12 m reduces to the plane-parallel run (none) for the same sky
// photons: the same rows, acceptance within 1e-5 relative, positions within 1e-3 m.
void test_spherical_limit_matches_none() {
  const auto photons = sky_photons(300);
  const auto none = records_of(propagate(sky_stream(run_keys(no_key), photons, 300),
                                         config_of(Curvature::none, kEarthRadius))
                                   .bytes);
  const auto spherical = records_of(
      propagate(sky_stream(run_keys(0), photons, 300), config_of(Curvature::spherical, 1e12))
          .bytes);
  for (const std::uint32_t telescope : {1U, 2U}) {
    const auto x_none = photon_column(none, telescope, Column::x);
    const auto x_sph = photon_column(spherical, telescope, Column::x);
    const auto y_none = photon_column(none, telescope, Column::y);
    const auto y_sph = photon_column(spherical, telescope, Column::y);
    const auto p_none = photon_column(none, telescope, Column::p);
    const auto p_sph = photon_column(spherical, telescope, Column::p);
    CHECK(!p_none.empty());
    CHECK(x_none.size() == x_sph.size() && y_none.size() == y_sph.size() &&
          p_none.size() == p_sph.size());
    if (x_none.size() != x_sph.size() || p_none.size() != p_sph.size() ||
        y_none.size() != y_sph.size()) {
      continue;
    }
    for (std::size_t i = 0; i < p_none.size(); ++i) {
      CHECK_NEAR(static_cast<double>(x_sph[i]), static_cast<double>(x_none[i]), 1e-3);
      CHECK_NEAR(static_cast<double>(y_sph[i]), static_cast<double>(y_none[i]), 1e-3);
      const auto flat = static_cast<double>(p_none[i]);
      CHECK_NEAR(static_cast<double>(p_sph[i]), flat, 1e-5 * std::abs(flat) + 1e-12);
    }
  }
}

// Emission altitude of a photon in the sky frame at lateral distance d: |p - C| - R with C = (0, 0,
// -(R + h)). For d^2 << R^2 the shift against the flat value a = h + z is d^2 / (2 (R + a)).
void test_sky_emission_altitude() {
  const double R = kEarthRadius;
  const std::vector<SkyPhoton> photon{{1500.0, 1500.0, 5000.0, 0.0, 0.0, 1.0, false}};
  const auto curved = records_of(
      propagate(sky_stream(run_keys(0), photon, 1), config_of(Curvature::spherical, R)).bytes);
  const auto flat = records_of(
      propagate(sky_stream(run_keys(no_key), photon, 1), config_of(Curvature::none, R)).bytes);
  const auto zem_curved = photon_column(curved, 1, Column::zem);
  const auto zem_flat = photon_column(flat, 1, Column::zem);
  CHECK(zem_curved.size() == 1 && zem_flat.size() == 1);
  if (zem_curved.size() != 1 || zem_flat.size() != 1) {
    return;
  }
  const double a = site_m + 5000.0;
  const double d2 = 1500.0 * 1500.0 * 2.0;
  const double exact = std::sqrt(d2 + (R + a) * (R + a)) - R;
  const double shift = static_cast<double>(zem_curved[0]) - static_cast<double>(zem_flat[0]);
  CHECK_NEAR(static_cast<double>(zem_flat[0]), a, 1e-2);
  CHECK_NEAR(static_cast<double>(zem_curved[0]), exact, 1e-2);
  CHECK(std::abs(shift - d2 / (2.0 * R)) <= 0.01 * d2 / (2.0 * R));
  std::fprintf(stderr, "sky zem shift %.6f m, d^2/(2R) %.6f m\n", shift, d2 / (2.0 * R));
}

// Curved modes: weighted acceptance per photon, thinned selection, and the output of a
// reversed row order do not depend on the chunk size or the order of the rows.
void test_chunk_and_order_independence() {
  const double R = 6.0e6;
  const auto photons = ground_photons(600);
  auto reversed = photons;
  std::reverse(reversed.begin(), reversed.end());

  for (const Curvature curvature : {Curvature::corsika, Curvature::spherical}) {
    const std::int64_t key = curvature == Curvature::corsika ? 1 : 0;
    const auto whole = records_of(
        propagate(ground_stream(run_keys(key), photons, 600), config_of(curvature, R)).bytes);
    const auto small = records_of(
        propagate(ground_stream(run_keys(key), photons, 7), config_of(curvature, R)).bytes);
    const auto rev = records_of(
        propagate(ground_stream(run_keys(key), reversed, 600), config_of(curvature, R)).bytes);
    const auto p_whole = value_by_x(whole, 1, Column::p);
    CHECK(p_whole.size() > 500);
    CHECK(p_whole == value_by_x(small, 1, Column::p));
    CHECK(p_whole == value_by_x(rev, 1, Column::p));

    auto thinned = config_of(curvature, R);
    thinned.mode = OutputMode::thinned;
    thinned.store_p = false;
    const auto thinned_whole =
        records_of(propagate(ground_stream(run_keys(key), photons, 600), thinned).bytes);
    const auto thinned_small =
        records_of(propagate(ground_stream(run_keys(key), photons, 7), thinned).bytes);
    CHECK(!sorted_x(thinned_whole, 1).empty());
    CHECK(sorted_x(thinned_whole, 1) == sorted_x(thinned_small, 1));
  }

  // Sky frame (spherical): per-photon geometry, chunk size independent.
  const auto sky = sky_photons(300);
  const auto sky_whole = records_of(
      propagate(sky_stream(run_keys(0), sky, 300), config_of(Curvature::spherical, R)).bytes);
  const auto sky_small = records_of(
      propagate(sky_stream(run_keys(0), sky, 7), config_of(Curvature::spherical, R)).bytes);
  CHECK(value_by_x(sky_whole, 1, Column::p) == value_by_x(sky_small, 1, Column::p));
  CHECK(value_by_x(sky_whole, 1, Column::zem) == value_by_x(sky_small, 1, Column::zem));
}

// Spherical refraction and air delay with R = 1e12 m reduce to the flat tables of mode none: the
// same rows, positions within 2 mm, directions within 1e-6, times within 1e-3 ns and acceptance
// within 1e-5 relative. Telescopes at the reference and 1 km east; emission 5 and 20 km above the
// plane at zenith angles 0 to 70 degrees.
void test_profile_limit_matches_none() {
  const std::vector<double> lateral{0.0, 1000.0};
  const std::vector<double> zero(2, 0.0);
  std::vector<SkyPhoton> photons;
  for (const double zenith : {0.0, 30.0, 60.0, 70.0}) {
    for (const double height : {5000.0, 20000.0}) {
      for (const double phi : {0.0, 0.7}) {
        for (const double x : {-300.0, 0.0, 300.0}) {
          photons.push_back(sky_at(x, 0.0, height, zenith, phi));
        }
      }
    }
  }
  const std::size_t n = photons.size();
  const auto none =
      records_of(propagate(sky_stream(lateral_run_keys(no_key, 0, lateral, zero), photons, n),
                           profile_config(Curvature::none, kEarthRadius))
                     .bytes);
  const auto spherical =
      records_of(propagate(sky_stream(lateral_run_keys(0, 0, lateral, zero), photons, n),
                           profile_config(Curvature::spherical, 1e12))
                     .bytes);
  for (const std::uint32_t telescope : {1U, 2U}) {
    const auto compare = [&](Column c, double abs_tol, double rel_tol) {
      const auto flat = photon_column(none, telescope, c);
      const auto curved = photon_column(spherical, telescope, c);
      // Photons of zero acceptance are removed by the weighted output in both runs alike.
      CHECK(!flat.empty());
      CHECK(curved.size() == flat.size());
      if (curved.size() != flat.size()) {
        return;
      }
      for (std::size_t i = 0; i < flat.size(); ++i) {
        const auto e = static_cast<double>(flat[i]);
        CHECK_NEAR(static_cast<double>(curved[i]), e, abs_tol + rel_tol * std::abs(e));
      }
    };
    compare(Column::x, 2e-3, 0.0);
    compare(Column::y, 2e-3, 0.0);
    compare(Column::cx, 1e-6, 0.0);
    compare(Column::cy, 1e-6, 0.0);
    compare(Column::t, 1e-3, 0.0);
    compare(Column::p, 1e-12, 1e-5);
  }
}

// Photons of the real-radius and off-grid checks (event 1): indices 0 (zenith 0) and 1 (zenith 80
// degrees) at 5 km, and 2 to 7 (zenith 60, 80, 85 degrees from 7 and 20 km), all on the axis of
// the telescope at the reference point.
std::vector<SkyPhoton> real_radius_photons() {
  std::vector<SkyPhoton> photons{sky_at(0.0, 0.0, 5000.0, 0.0), sky_at(0.0, 0.0, 5000.0, 80.0)};
  for (const double zenith : {60.0, 80.0, 85.0}) {
    for (const double height : {7000.0, 20000.0}) {
      photons.push_back(sky_at(0.0, 0.0, height, zenith));
    }
  }
  return photons;
}

// Output of one propagation with the real Earth radius and the profile (one table cache for all
// checks): run A (refraction not applied, telescopes 1 to 4 at z = 0, 0, 3.7, 9.9 m) with events
// 1 to 6, run B (refraction applied) with event 7. Geometry only (no extinction), so the rows of
// the output follow the input order, apart from the dropped photons.
struct Batch {
  obdeect::propagation::PropagationSummary summary{};
  std::vector<std::vector<obdeect::stream::Record>> events;
  std::vector<obdeect::stream::Record> records;
};

Batch real_radius_batch() {
  const std::vector<double> zero(4, 0.0);
  const auto run_a = lateral_run_keys(0, 0, {0.0, 1000.0, 0.0, 0.0}, zero, {0.0, 0.0, 3.7, 9.9});
  const auto run_b = lateral_run_keys(0, 1, {0.0, 1000.0, 0.0, 0.0}, zero, zero);
  const auto photons = real_radius_photons();
  // Chunk and order checks: explicit wavelengths, so reordering does not change the draws.
  auto sorted = sky_photons(300);
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    sorted[i].wavelength = 300.0 + static_cast<double>(i % 400);
  }
  auto reversed = sorted;
  std::reverse(reversed.begin(), reversed.end());
  const std::vector<SkyPhoton> dropped{sky_at(0.0, 0.0, 5000.0, 0.0),
                                       sky_at(0.0, 0.0, 17000.0, 89.0)};
  const std::string bytes = sky_blocks({
      {run_a, 1, photons, photons.size(), 4},             // real radius and off-grid telescopes
      {run_a, 2, sorted, 300, 2},                         // whole chunk
      {run_a, 3, sorted, 7, 2},                           // chunks of 7 rows
      {run_a, 4, reversed, 300, 2},                       // reversed order
      {run_a, 5, sorted, 300, 2},                         // the same input again (repeat)
      {run_a, 6, dropped, 2, 2},                          // one photon above the telescope level
      {run_b, 7, {sky_at(0.0, 0.0, 4000.0, 60.0)}, 1, 1}, // refraction applied in the input
  });
  auto config = profile_config(Curvature::spherical, kEarthRadius);
  config.geometry_only = true;
  config.transmission.reset();
  config.store_p = false;
  const auto result = propagate(bytes, config);
  Batch batch;
  batch.summary = result.summary;
  batch.records = records_of(result.bytes);
  batch.events = split_events(batch.records);
  CHECK(batch.events.size() == 8);
  return batch;
}

// Real Earth radius with a profile (refraction_applied = 0): on the telescope axis, zenith 0 equals
// the flat result; at zenith 80 degrees the arrival point is the straight-line position plus the
// shift of curved_ray_integrals, and the direction is scaled by (1 + nu_e) / (1 + nu_t).
void check_real_radius(const Batch &batch) {
  const double R = kEarthRadius;
  const auto profile = ecmwf_profile();
  const std::vector<double> zero(2, 0.0);
  const auto photons = real_radius_photons();
  auto flat_config = profile_config(Curvature::none, R);
  flat_config.geometry_only = true;
  flat_config.transmission.reset();
  flat_config.store_p = false;
  const std::string flat_bytes = sky_blocks(
      {{lateral_run_keys(no_key, 0, {0.0, 1000.0}, zero), 1, {photons[0], photons[1]}, 2, 1}});
  const auto none = records_of(propagate(flat_bytes, flat_config).bytes);
  const auto &event = batch.events[1];
  const auto x_flat = photon_column(none, 1, Column::x);
  const auto x_curved = photon_column(event, 1, Column::x);
  CHECK(x_flat.size() == 2 && x_curved.size() == photons.size());
  if (x_flat.size() != 2 || x_curved.size() != photons.size()) {
    return;
  }

  // Zenith 0: no bending, the same air delay.
  CHECK_NEAR(static_cast<double>(x_curved[0]), static_cast<double>(x_flat[0]), 1e-3);
  CHECK_NEAR(static_cast<double>(photon_column(event, 1, Column::y)[0]),
             static_cast<double>(photon_column(none, 1, Column::y)[0]), 1e-3);
  CHECK_NEAR(static_cast<double>(photon_column(event, 1, Column::t)[0]),
             static_cast<double>(photon_column(none, 1, Column::t)[0]), 1e-3);

  // Zenith 80 degrees: the chunk stores cx as float32, the reference uses the same value.
  const double cx = static_cast<double>(static_cast<float>(photons[1].cx));
  const double cz = std::sqrt(1.0 - cx * cx);
  const double straight = 5000.0 / cz * cx; // telescope 1 at the reference point
  const double emission_m = site_m + 5000.0;
  const double r_t = R + site_m;
  const double b =
      obdeect::propagation::impact_parameter_m(0.0, 0.0, 5000.0, site_m, R, cx, 0.0, cz);
  const double theta = std::atan2(b, std::sqrt((r_t - b) * (r_t + b)));
  const auto direct =
      obdeect::propagation::curved_ray_integrals(profile, R, site_m, theta, emission_m);
  CHECK_NEAR(static_cast<double>(x_curved[1]), straight + direct.shift_m, 1e-3);
  const double shift_flat = static_cast<double>(x_flat[1]) - straight;
  const double shift_curved = static_cast<double>(x_curved[1]) - straight;
  std::fprintf(stderr, "zenith 80 deg: shift flat %.4f m, spherical %.4f m\n", shift_flat,
               shift_curved);
  CHECK(std::abs(shift_curved - shift_flat) > 0.01 * std::abs(shift_flat));

  const double nu_e = profile.index_minus_one_cm(100.0 * emission_m);
  const double nu_t = profile.index_minus_one_cm(100.0 * site_m);
  const auto cx_out = photon_column(event, 1, Column::cx);
  CHECK(cx_out.size() == photons.size());
  if (cx_out.size() == photons.size()) {
    CHECK_NEAR(static_cast<double>(cx_out[1]) / cx, (1.0 + nu_e) / (1.0 + nu_t), 1e-6);
  }
  const auto zem = photon_column(event, 1, Column::zem);
  CHECK(zem.size() == photons.size());
  if (zem.size() == photons.size()) {
    CHECK_NEAR(static_cast<double>(zem[1]), emission_m, 1e-2);
  }
}

// refraction_applied = 1 in the input (run B): no bending, only the air delay of the straight path.
void check_refraction_applied(const Batch &batch) {
  const double R = kEarthRadius;
  const auto profile = ecmwf_profile();
  const double cx = static_cast<double>(static_cast<float>(sky_at(0.0, 0.0, 4000.0, 60.0).cx));
  const double cz = std::sqrt(1.0 - cx * cx);
  const double path = 4000.0 / cz;
  const double b =
      obdeect::propagation::impact_parameter_m(0.0, 0.0, 4000.0, site_m, R, cx, 0.0, cz);
  const double r_t = R + site_m;
  const double theta = std::atan2(b, std::sqrt((r_t - b) * (r_t + b)));
  const auto direct =
      obdeect::propagation::curved_ray_integrals(profile, R, site_m, theta, site_m + 4000.0);
  const auto &event = batch.events[7];
  const auto x = photon_column(event, 1, Column::x);
  const auto t = photon_column(event, 1, Column::t);
  const auto cx_out = photon_column(event, 1, Column::cx);
  CHECK(x.size() == 1 && t.size() == 1 && cx_out.size() == 1);
  if (x.size() == 1 && t.size() == 1 && cx_out.size() == 1) {
    CHECK_NEAR(static_cast<double>(x[0]), path * cx, 1e-3);
    CHECK_NEAR(static_cast<double>(cx_out[0]), cx, 1e-7);
    CHECK_NEAR(static_cast<double>(t[0]),
               path / obdeect::propagation::speed_of_light_m_per_ns + direct.excess_straight_ns,
               1e-3);
    CHECK(direct.excess_straight_ns > 0.0);
  }
  // The output run record of run B says that refraction is applied.
  const obdeect::stream::Metadata *last_run = nullptr;
  for (const auto &record : batch.records) {
    if (record.type == RecordType::run) {
      last_run = &record.metadata;
    }
  }
  CHECK(last_run != nullptr);
  if (last_run != nullptr) {
    CHECK(std::get<std::int64_t>(last_run->at("refraction_applied")) == 1);
  }
}

// A photon whose line does not reach the telescope level (b >= r_t) is dropped with w = 0: zenith
// 89 degrees from 19 km altitude passes above the telescope at 2 km. Event 6 of the batch, and the
// total of the batch: the same drops as the upward photons of events 2 to 5 at two telescopes.
// Without a profile (no tables) the same photon is dropped in a separate run.
void check_dropped(const Batch &batch) {
  const auto &event = batch.events[6];
  for (const std::uint32_t telescope : {1U, 2U}) {
    const auto zem = photon_column(event, telescope, Column::zem);
    CHECK(zem.size() == 1);
    if (zem.size() == 1) {
      CHECK_NEAR(static_cast<double>(zem[0]), site_m + 5000.0, 1e-2);
    }
  }
  std::uint64_t upward = 0;
  for (const auto &photon : sky_photons(300)) {
    upward += photon.upward ? 1U : 0U;
  }
  CHECK(batch.summary.dropped == 2 + 4 * 2 * upward);

  const std::vector<double> lateral{0.0, 1000.0};
  const std::vector<double> zero(2, 0.0);
  const std::vector<SkyPhoton> photons{sky_at(0.0, 0.0, 5000.0, 0.0),
                                       sky_at(0.0, 0.0, 17000.0, 89.0)};
  const auto result =
      propagate(sky_blocks({{lateral_run_keys(0, 0, lateral, zero), 1, photons, 2, 2}}),
                config_of(Curvature::spherical, kEarthRadius));
  // Two telescopes, one dropped photon each; the weighted output keeps the other photon.
  CHECK(result.summary.dropped == 2);
  CHECK(result.summary.photons_in == 4);
  CHECK(result.summary.photons_out == 2);
  CHECK(photon_column(records_of(result.bytes), 1, Column::x).size() == 1);
}

// Spherical refraction with a profile: the same photons in whole chunks, chunks of 7 rows, reversed
// order and again (repeat) give the same output per photon; the upward photons are dropped.
void check_chunk_order(const Batch &batch) {
  std::uint64_t upward = 0;
  for (const auto &photon : sky_photons(300)) {
    upward += photon.upward ? 1U : 0U;
  }
  const auto whole = value_by_x(batch.events[2], 1, Column::t);
  CHECK(whole.size() == 300 - upward);
  CHECK(!whole.empty());
  for (const std::size_t event : {3U, 4U, 5U}) {
    for (const Column c : {Column::t, Column::zem, Column::cx, Column::cy}) {
      CHECK(value_by_x(batch.events[2], 1, c) == value_by_x(batch.events[event], 1, c));
    }
  }
}

// SkyPlane contract of sky_to_ground in the spherical geometry: a profile needs the two curved
// tables bracketing the telescope altitude, flat tables are refused, and the tables must match the
// plane's altitude, weight and radius.
void test_sky_plane_validation() {
  const double R = kEarthRadius;
  const auto profile = ecmwf_profile();
  const obdeect::propagation::RefractionTables flat(profile, site_m);
  const obdeect::propagation::CurvedRefractionTables lower(profile, R, site_m, false,
                                                           build_threads());
  const obdeect::propagation::CurvedRefractionTables upper(profile, R, site_m + 10.0, false,
                                                           build_threads());
  const auto sky_one = [] {
    return sky_chunk(1, 0, {SkyPhoton{0.0, 0.0, 5000.0, 0.0, 0.0, 1.0, false, 0.0}}, 1);
  };
  const auto plane = [&](bool with_profile) {
    obdeect::propagation::SkyPlane p;
    p.observation_m = site_m;
    p.earth_radius_m = R;
    p.telescope_altitude_m = site_m;
    if (with_profile) {
      p.profile = &profile;
    }
    return p;
  };
  const auto convert = [&](const obdeect::propagation::SkyPlane &p) {
    auto chunk = sky_one();
    (void)obdeect::propagation::sky_to_ground(chunk, p);
  };

  CHECK(throws_with([&] { convert(plane(true)); },
                    "lower curved table")); // profile without curved tables
  auto with_tables = plane(true);
  with_tables.curved_tables = &lower;
  with_tables.curved_tables_upper = &upper;
  with_tables.curved_bending = true;
  CHECK(accepts([&] { convert(with_tables); }));
  // The lower table alone is valid with weight 0 (telescope on the lower level, or no upper level
  // inside the profile) and invalid with a weight above 0.
  auto one_table = with_tables;
  one_table.curved_tables_upper = nullptr;
  CHECK(accepts([&] { convert(one_table); }));
  auto one_table_weighted = one_table;
  one_table_weighted.telescope_altitude_m = site_m + 3.0;
  one_table_weighted.curved_weight = 0.3;
  CHECK(throws_with([&] { convert(one_table_weighted); }, "bracket"));
  auto upper_only = with_tables;
  upper_only.curved_tables = nullptr;
  CHECK(throws_with([&] { convert(upper_only); }, "lower curved table"));
  // A photon emitted above the telescope but below the upper table level is kept: the lower table
  // alone serves it (the upper lookup starts at its own level), with zero or positive weight.
  const auto between_levels = [&](const obdeect::propagation::SkyPlane &p, double z_m) {
    auto chunk = sky_chunk(1, 0, {sky_at(0.0, 0.0, z_m, 30.0)}, 1);
    const auto lost = obdeect::propagation::sky_to_ground(chunk, p);
    CHECK(lost.total() == 0);
    CHECK(chunk[Column::w][0] > 0.0F);
    return chunk;
  };
  // Where no table is referenced to the telescope level the integrals are computed directly at the
  // telescope altitude: the output equals the straight line plus curved_ray_integrals there.
  const auto expect_direct = [&](const obdeect::propagation::SkyPlane &p, double z_m) {
    auto input = sky_chunk(1, 0, {sky_at(0.0, 0.0, z_m, 30.0)}, 1);
    const double cx = input[Column::cx][0];
    const double cy = input[Column::cy][0];
    const double cz = std::sqrt(1.0 - cx * cx - cy * cy);
    const double path = (z_m - p.telescope_z_m) / cz;
    const double b = obdeect::propagation::impact_parameter_m(0.0, 0.0, z_m, site_m, R, cx, cy, cz);
    const auto out = between_levels(p, z_m);
    const double theta = std::asin(b / (R + p.telescope_altitude_m));
    const auto direct = obdeect::propagation::curved_ray_integrals(
        profile, R, p.telescope_altitude_m, theta, static_cast<double>(out[Column::zem][0]));
    const double s = std::sqrt(cx * cx + cy * cy);
    CHECK(std::abs(out[Column::x][0] - (path * cx + direct.shift_m * cx / s)) < 1e-4);
    CHECK(std::abs(out[Column::t][0] - (path / obdeect::propagation::speed_of_light_m_per_ns +
                                        direct.excess_bent_ns)) < 2e-4);
  };
  CHECK(std::isfinite(between_levels(with_tables, 5.0)[Column::t][0])); // on the level, 5 m above
  auto off_level = with_tables;
  off_level.telescope_z_m = 3.0;
  off_level.telescope_altitude_m = site_m + 3.0;
  off_level.curved_weight = 0.3;
  expect_direct(off_level, 8.0); // emitted 5 m above the telescope, below site + 10
  auto flat_tables = plane(true);
  flat_tables.tables = &flat;
  CHECK(throws_with([&] { convert(flat_tables); }, "curved geometry takes curved tables"));
  auto wrong_altitude = with_tables;
  wrong_altitude.telescope_altitude_m = site_m + 100.0; // outside [lower, upper)
  CHECK(throws_with([&] { convert(wrong_altitude); }, "bracket"));
  auto wrong_weight = with_tables;
  wrong_weight.telescope_altitude_m = site_m + 3.0; // weight must be 0.3
  CHECK(throws_with([&] { convert(wrong_weight); }, "bracket"));
  auto no_profile = plane(false);
  no_profile.curved_tables = &lower;
  no_profile.curved_tables_upper = &upper;
  CHECK(throws_with([&] { convert(no_profile); }, "needs a profile"));
  auto flat_radius = plane(true);
  flat_radius.earth_radius_m = 0.0;
  flat_radius.curved_tables = &lower;
  flat_radius.curved_tables_upper = &upper;
  CHECK(throws_with([&] { convert(flat_radius); }, "Earth radius"));
}

// Telescope altitudes that are not on the 10 m grid (observation level 2000 m, telescope z = 3.7
// and 9.9 m, event 1 of the batch). Photons at zenith 60, 80 and 85 degrees from 7 and 20 km, on
// the telescope axis: x, y and t against curved_ray_integrals at the exact altitude plus the
// straight-line values. The time tolerance is 0.01 ns plus the float32 half ulp of the output t at
// that value (the output t column is float32, not changed here) plus 5e-6 of the bent excess, the
// design accuracy of the linear interpolation in the 10 m spacing (measured 2e-6 to 3e-6 at 85
// degrees); the rounding of the output and the interpolation error add up. The position tolerances
// are 1 cm at 60 and 80 degrees, 5 cm at 85.
void check_off_grid(const Batch &batch) {
  const double R = kEarthRadius;
  const auto profile = ecmwf_profile();
  const double c = obdeect::propagation::speed_of_light_m_per_ns;
  const auto photons = real_radius_photons();
  const auto &event = batch.events[1];
  for (const std::uint32_t telescope : {3U, 4U}) {
    const double z = telescope == 3U ? 3.7 : 9.9;
    const auto x_out = photon_column(event, telescope, Column::x);
    const auto y_out = photon_column(event, telescope, Column::y);
    const auto t_out = photon_column(event, telescope, Column::t);
    CHECK(x_out.size() == photons.size() && y_out.size() == x_out.size() &&
          t_out.size() == x_out.size());
    if (x_out.size() != photons.size() || y_out.size() != photons.size() ||
        t_out.size() != photons.size()) {
      continue;
    }
    const double h_t = site_m + z; // exact telescope altitude
    for (std::size_t i = 2; i < photons.size(); ++i) {
      const SkyPhoton &photon = photons[i];
      const double zenith = i < 4 ? 60.0 : (i < 6 ? 80.0 : 85.0);
      const double cx = static_cast<double>(static_cast<float>(photon.cx));
      const double cz = std::sqrt(1.0 - cx * cx);
      const double path = (photon.z - z) / cz;
      const double b =
          obdeect::propagation::impact_parameter_m(0.0, 0.0, photon.z, site_m, R, cx, 0.0, cz);
      const double r_t = R + h_t;
      const double theta = std::atan2(b, std::sqrt((r_t - b) * (r_t + b)));
      const auto direct =
          obdeect::propagation::curved_ray_integrals(profile, R, h_t, theta, site_m + photon.z);
      const double x_expected = path * cx + direct.shift_m;
      const double t_expected = path / c + direct.excess_bent_ns;
      const float t_float = static_cast<float>(t_expected);
      const double half_ulp =
          0.5 *
          (static_cast<double>(std::nextafter(t_float, std::numeric_limits<float>::infinity())) -
           static_cast<double>(t_float));
      const double pos_tol = zenith == 85.0 ? 5e-2 : 1e-2;
      const double time_tol = 0.01 + half_ulp + 5e-6 * std::abs(direct.excess_bent_ns);
      const double ex = std::abs(static_cast<double>(x_out[i]) - x_expected);
      const double ey = std::abs(static_cast<double>(y_out[i]));
      const double et = std::abs(static_cast<double>(t_out[i]) - t_expected);
      std::fprintf(stderr,
                   "z %.1f m, zenith %2.0f, H %6.0f m: |dx| %.2e m (tol %.0e), |dy| %.2e m, "
                   "|dt| %.2e ns (tol %.2e)\n",
                   z, zenith, photon.z, ex, pos_tol, ey, et, time_tol);
      CHECK(ex <= pos_tol);
      CHECK(ey <= pos_tol);
      CHECK(et <= time_tol);
    }
  }
}

} // namespace

// The telescope level lies in the last 10 m below the profile top, so no upper table level fits
// inside the profile: the propagator uses the lower table alone and keeps the photons.
void test_profile_top_boundary() {
  const double top = ecmwf_profile().top_m();
  const double level = top - 5.0;
  auto run = lateral_run_keys(0, 0, {0.0}, {0.0});
  run["observation_level_m"] = level;
  const std::vector<SkyPhoton> photons{sky_at(0.0, 0.0, 3.0, 30.0), sky_at(0.0, 0.0, 4.0, 60.0)};
  auto config = profile_config(Curvature::spherical, kEarthRadius);
  config.geometry_only = true;
  config.transmission.reset();
  config.store_p = false;
  config.site_altitude_m = level;
  const auto result = propagate(sky_blocks({{run, 1, photons, photons.size(), 1}}), config);
  CHECK(result.summary.dropped == 0);
  CHECK(result.summary.photons_out == photons.size());
}

int main() {
  return run_checks([] {
    test_rejections();
    test_status_and_run_key();
    test_corsika_extinction();
    test_spherical_limit_matches_none();
    test_sky_emission_altitude();
    test_chunk_and_order_independence();
    test_profile_limit_matches_none();
    // One propagation with the real radius and a profile; its tables serve all checks below.
    const Batch batch = real_radius_batch();
    check_real_radius(batch);
    check_off_grid(batch);
    check_refraction_applied(batch);
    check_dropped(batch);
    check_chunk_order(batch);
    test_sky_plane_validation();
    test_profile_top_boundary();
  });
}
