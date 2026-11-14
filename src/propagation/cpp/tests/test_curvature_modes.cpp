// Curvature modes of the atmosphere stage (propagator): input checks and rejections, efficiency
// status keys, the CORSIKA CURVED extinction, the spherical sky geometry (limit R -> infinity and
// the Earth-radius correction), and independence of chunk size and row order in the curved modes.
#include "check.hpp"

#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/propagator.hpp"
#include "obdeect/propagation/scattering.hpp"
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
#include <variant>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

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

std::string data_file(const std::string &name) {
  return std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/" + name;
}

obdeect::propagation::TransmissionTable real_table() {
  return obdeect::propagation::TransmissionTable::from_file(
      data_file("atm_trans_2000_1_10_0_0_2000.dat"), site_m);
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

obdeect::stream::Metadata event_keys() {
  obdeect::stream::Metadata m;
  m["event_number"] = std::int64_t{7};
  return m;
}

obdeect::stream::Metadata telescope_keys(std::uint32_t id, Frame frame) {
  obdeect::stream::Metadata m;
  m["telescope_id"] = static_cast<std::int64_t>(id);
  m["frame"] = static_cast<std::int64_t>(frame);
  return m;
}

obdeect::stream::Metadata telescope_end_keys(std::uint32_t id) {
  obdeect::stream::Metadata m;
  m["telescope_id"] = static_cast<std::int64_t>(id);
  m["photons"] = 0.0;
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
  const std::uint32_t mask =
      obdeect::stream::column_bit(Column::x) | obdeect::stream::column_bit(Column::y) |
      obdeect::stream::column_bit(Column::z) | obdeect::stream::column_bit(Column::cx) |
      obdeect::stream::column_bit(Column::cy) | obdeect::stream::column_bit(Column::cz) |
      obdeect::stream::column_bit(Column::t) | obdeect::stream::column_bit(Column::w);
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

std::vector<obdeect::stream::Record> records_of(const std::string &bytes) {
  std::istringstream in(bytes);
  obdeect::stream::Reader reader(in);
  std::vector<obdeect::stream::Record> records;
  obdeect::stream::Record record;
  while (reader.next(record)) {
    records.push_back(record);
  }
  CHECK(reader.ended());
  return records;
}

// Values of column c over all photon records of one telescope.
std::vector<float> photon_column(const std::vector<obdeect::stream::Record> &records,
                                 std::uint32_t telescope, Column c) {
  std::vector<float> values;
  for (const auto &record : records) {
    if (record.type == RecordType::photons && record.chunk.telescope_id == telescope &&
        record.chunk.has(c)) {
      values.insert(values.end(), record.chunk[c].begin(), record.chunk[c].end());
    }
  }
  return values;
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

bool throws_invalid(const std::function<void()> &fn) {
  try {
    fn();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
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
  CHECK(throws_invalid([&] { propagate(ground_in(1), config_of(Curvature::none, R)); }));
  CHECK(!throws_invalid([&] { propagate(ground_in(no_key), config_of(Curvature::none, R)); }));
  CHECK(!throws_invalid([&] { propagate(ground_in(0), config_of(Curvature::none, R)); }));

  // corsika: the key must be present and 1.
  CHECK(throws_invalid([&] { propagate(ground_in(no_key), config_of(Curvature::corsika, R)); }));
  CHECK(throws_invalid([&] { propagate(ground_in(0), config_of(Curvature::corsika, R)); }));
  CHECK(!throws_invalid([&] { propagate(ground_in(1), config_of(Curvature::corsika, R)); }));
  // corsika: sky-frame telescopes have no CORSIKA CURVED geometry.
  CHECK(throws_invalid([&] { propagate(sky_in(1), config_of(Curvature::corsika, R)); }));

  // spherical: the key must be present and 0.
  CHECK(throws_invalid([&] { propagate(ground_in(no_key), config_of(Curvature::spherical, R)); }));
  CHECK(throws_invalid([&] { propagate(ground_in(1), config_of(Curvature::spherical, R)); }));
  CHECK(!throws_invalid([&] { propagate(sky_in(0), config_of(Curvature::spherical, R)); }));

  // Scattering and a profile are not combined with the curved modes (checked at construction).
  for (const Curvature curvature : {Curvature::corsika, Curvature::spherical}) {
    auto scattering = config_of(curvature, R);
    scattering.scattering.enabled = true;
    scattering.scattering.points_per_photon = 8;
    scattering.scattering_table.emplace(
        obdeect::propagation::ScatteringTable::from_file(data_file("scattering_exponential.dat")));
    CHECK(throws_invalid([&] {
      obdeect::propagation::Propagator propagator(scattering);
      (void)propagator;
    }));
  }
  auto profile = config_of(Curvature::spherical, R);
  profile.profile.emplace(obdeect::propagation::AtmosphereProfile::from_file(
      data_file("atmprof_ecmwf_north_winter_fixed.dat")));
  profile.profile_source = "fixture";
  CHECK(throws_invalid([&] {
    obdeect::propagation::Propagator propagator(profile);
    (void)propagator;
  }));

  // Earth radius: finite and above 1e6 m, in every mode.
  for (const Curvature curvature : {Curvature::none, Curvature::corsika, Curvature::spherical}) {
    for (const double radius : {0.0, -1.0, 1e6, std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::infinity()}) {
      CHECK(throws_invalid([&] {
        obdeect::propagation::Propagator propagator(config_of(curvature, radius));
        (void)propagator;
      }));
    }
  }

  // Mode names.
  CHECK(obdeect::propagation::parse_curvature("none") == Curvature::none);
  CHECK(obdeect::propagation::parse_curvature("corsika") == Curvature::corsika);
  CHECK(obdeect::propagation::parse_curvature("spherical") == Curvature::spherical);
  for (const char *name : {"", "Corsika", "curved", "flat"}) {
    CHECK(throws_invalid([&] { obdeect::propagation::parse_curvature(name); }));
  }

  // Geometry only: the curved geometry is allowed (no extinction).
  auto geometry_spherical = config_of(Curvature::spherical, R);
  geometry_spherical.geometry_only = true;
  geometry_spherical.transmission.reset();
  geometry_spherical.store_p = false;
  CHECK(!throws_invalid([&] { propagate(sky_in(0), geometry_spherical); }));
  auto geometry_corsika = config_of(Curvature::corsika, R);
  geometry_corsika.geometry_only = true;
  geometry_corsika.transmission.reset();
  geometry_corsika.store_p = false;
  CHECK(!throws_invalid([&] { propagate(ground_in(1), geometry_corsika); }));
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

// Emission altitude of a sky photon at lateral distance d: |p - C| - R with C = (0, 0, -(R + h)).
// For d^2 << R^2 the shift against the flat value a = h + z is d^2 / (2 (R + a)).
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

} // namespace

int main() {
  return run_checks([] {
    test_rejections();
    test_status_and_run_key();
    test_corsika_extinction();
    test_spherical_limit_matches_none();
    test_sky_emission_altitude();
    test_chunk_and_order_independence();
  });
}
