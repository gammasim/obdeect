#include "check.hpp"

#include "dlufd/profile.hpp"
#include "dlufd/propagator.hpp"
#include "dlufd/refraction.hpp"
#include "dlufd/scattering.hpp"
#include "dlufd/timing.hpp"
#include "dlufd/transmission.hpp"

#include "gstrickd/photons.hpp"
#include "gstrickd/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#ifndef DLUFD_TEST_DATA_DIR
#error "DLUFD_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using gstrickd::Column;
using gstrickd::Frame;
using gstrickd::RecordType;

constexpr double site_m = 2000.0;

std::int64_t i64(std::int64_t value) { return value; }

// One photon of a synthetic chunk. Sky photons use z, cz and t; G photons use zem.
struct Row {
  double x = 0.0, y = 0.0, z = 0.0, cx = 0.0, cy = 0.0, cz = 0.0, t = 0.0, zem = 0.0, w = 1.0;
};

// Transmission without absorption (optical depth 0) for 200..999 nm, at the site altitude.
dlufd::TransmissionTable zero_table() {
  std::string text = "# H2= 2.000, H1=    3.000     5.000\n";
  for (int nm = 200; nm < 1000; ++nm) {
    text += std::to_string(nm) + " 0 0\n";
  }
  return dlufd::TransmissionTable(text, site_m);
}

dlufd::TransmissionTable real_table() {
  return dlufd::TransmissionTable::from_file(
      std::string(DLUFD_TEST_DATA_DIR) + "/atm_trans_2000_1_10_0_0_2000.dat", site_m);
}

dlufd::AtmosphereProfile ecmwf_profile() {
  return dlufd::AtmosphereProfile::from_file(std::string(DLUFD_TEST_DATA_DIR) +
                                             "/atmprof_ecmwf_north_winter_fixed.dat");
}

// Run keys: two telescopes at x = 100 m and x = -50 m; the observation level at 2000 m.
// refraction < 0: the key is absent.
gstrickd::Metadata run_keys(int refraction) {
  gstrickd::Metadata m;
  m["run_number"] = i64(1);
  m["observation_level_m"] = site_m;
  m["telescope_x_m"] = std::vector<double>{100.0, -50.0};
  m["telescope_y_m"] = std::vector<double>{0.0, 300.0};
  m["telescope_z_m"] = std::vector<double>{0.0, 0.0};
  m["telescope_r_m"] = std::vector<double>{12.0, 12.0};
  if (refraction >= 0) {
    m["refraction_applied"] = i64(refraction);
  }
  m["wavelength_min_nm"] = 250.0;
  m["wavelength_max_nm"] = 700.0;
  return m;
}

gstrickd::Metadata event_keys() {
  gstrickd::Metadata m;
  m["event_number"] = i64(7);
  return m;
}

gstrickd::Metadata telescope_keys(std::uint32_t id, Frame frame) {
  gstrickd::Metadata m;
  m["telescope_id"] = i64(id);
  m["frame"] = i64(static_cast<std::int64_t>(frame));
  return m;
}

gstrickd::Metadata telescope_end_keys(std::uint32_t id) {
  gstrickd::Metadata m;
  m["telescope_id"] = i64(id);
  m["photons"] = 0.0;
  return m;
}

// Photons of one chunk: sky (x y z cx cy cz t w) or G (x y cx cy t zem w).
gstrickd::PhotonChunk make_chunk(Frame frame, std::uint32_t telescope, std::uint64_t first,
                                 std::size_t n, const std::function<Row(std::size_t)> &row_of) {
  gstrickd::PhotonChunk chunk;
  chunk.telescope_id = telescope;
  chunk.frame = frame;
  chunk.first_index = first;
  const bool sky = frame == Frame::sky;
  std::uint32_t mask = gstrickd::column_bit(Column::x) | gstrickd::column_bit(Column::y) |
                       gstrickd::column_bit(Column::cx) | gstrickd::column_bit(Column::cy) |
                       gstrickd::column_bit(Column::t) | gstrickd::column_bit(Column::w);
  mask |= sky ? gstrickd::column_bit(Column::z) | gstrickd::column_bit(Column::cz)
              : gstrickd::column_bit(Column::zem);
  chunk.resize(n, mask);
  for (std::size_t r = 0; r < n; ++r) {
    const Row row = row_of(static_cast<std::size_t>(first) + r);
    chunk[Column::x][r] = static_cast<float>(row.x);
    chunk[Column::y][r] = static_cast<float>(row.y);
    chunk[Column::cx][r] = static_cast<float>(row.cx);
    chunk[Column::cy][r] = static_cast<float>(row.cy);
    chunk[Column::t][r] = static_cast<float>(row.t);
    chunk[Column::w][r] = static_cast<float>(row.w);
    if (sky) {
      chunk[Column::z][r] = static_cast<float>(row.z);
      chunk[Column::cz][r] = static_cast<float>(row.cz);
    } else {
      chunk[Column::zem][r] = static_cast<float>(row.zem);
    }
  }
  return chunk;
}

// G-frame photons with a spread of positions, directions, emission altitudes and weights.
Row ground_row(std::size_t i) {
  Row row;
  row.x = 300.0 * std::sin(0.37 * static_cast<double>(i));
  row.y = 300.0 * std::cos(0.23 * static_cast<double>(i));
  row.cx = 0.4 * std::sin(0.001 * static_cast<double>(i));
  row.cy = 0.4 * std::cos(0.0013 * static_cast<double>(i));
  row.t = static_cast<double>(i % 13);
  row.zem = 2000.0 + 1000.0 * static_cast<double>(i % 7);
  row.w = 1.0 + static_cast<double>(i % 4);
  return row;
}

// Sky photons, emitted 3 to 23 km above the observation level, some moving upwards.
Row sky_row(std::size_t i) {
  Row row;
  row.x = 500.0 * std::sin(0.1 * static_cast<double>(i));
  row.y = 500.0 * std::cos(0.07 * static_cast<double>(i));
  row.z = 3000.0 + 5000.0 * static_cast<double>(i % 5);
  row.cx = 0.3 * std::sin(0.01 * static_cast<double>(i));
  row.cy = 0.3 * std::cos(0.013 * static_cast<double>(i));
  const double magnitude = std::sqrt(1.0 - row.cx * row.cx - row.cy * row.cy);
  row.cz = i % 17 == 0 ? magnitude : -magnitude;
  row.t = static_cast<double>(i % 11);
  row.w = 1.0 + static_cast<double>(i % 3);
  return row;
}

std::string make_stream(const std::function<void(gstrickd::Writer &)> &body) {
  std::ostringstream out;
  {
    gstrickd::Writer writer(out, "test");
    body(writer);
    writer.finish();
  }
  return out.str();
}

std::vector<gstrickd::Record> records_of(const std::string &bytes) {
  std::istringstream in(bytes);
  gstrickd::Reader reader(in);
  std::vector<gstrickd::Record> records;
  gstrickd::Record record;
  while (reader.next(record)) {
    records.push_back(record);
  }
  CHECK(reader.ended());
  return records;
}

// Concatenated values of column c over all photon records of one telescope.
std::vector<float> photon_column(const std::vector<gstrickd::Record> &records,
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

std::size_t photon_rows(const std::vector<gstrickd::Record> &records, std::uint32_t telescope) {
  std::size_t rows = 0;
  for (const auto &record : records) {
    if (record.type == RecordType::photons && record.chunk.telescope_id == telescope) {
      rows += record.chunk.rows();
    }
  }
  return rows;
}

struct Result {
  std::string bytes;
  dlufd::PropagationSummary summary;
};

Result propagate(const std::string &input, const dlufd::PropagatorConfig &config) {
  std::istringstream in(input);
  std::ostringstream out;
  Result result;
  result.summary = dlufd::propagate_stream(in, out, config);
  result.bytes = out.str();
  return result;
}

dlufd::PropagatorConfig config_for(const dlufd::TransmissionTable &table) {
  dlufd::PropagatorConfig config;
  config.transmission.emplace(table);
  config.transmission_source = "test table";
  config.site_altitude_m = site_m;
  config.seed = 20261107;
  return config;
}

dlufd::EnvelopeFactor constant_factor(const std::string &name, double value) {
  return dlufd::EnvelopeFactor{name, {200.0, 800.0}, {value, value}, name + ".dat"};
}

// Stream with telescopes 1 and 2 of one event: G-frame photons, n per telescope.
std::string ground_stream(std::size_t n, std::size_t chunk_rows, int refraction = 0) {
  return make_stream([&](gstrickd::Writer &writer) {
    writer.write_metadata(RecordType::run, run_keys(refraction));
    writer.write_metadata(RecordType::event, event_keys());
    for (const std::uint32_t telescope : {1U, 2U}) {
      writer.write_metadata(RecordType::telescope, telescope_keys(telescope, Frame::ground));
      for (std::size_t first = 0; first < n; first += chunk_rows) {
        const std::size_t rows = std::min(chunk_rows, n - first);
        writer.write_chunk(make_chunk(Frame::ground, telescope, first, rows, ground_row));
      }
      writer.write_metadata(RecordType::telescope_end, telescope_end_keys(telescope));
    }
    writer.write_metadata(RecordType::event_end, {});
    writer.write_metadata(RecordType::run_end, {});
  });
}

// Ledger record that follows the run record (first one), or nullptr.
const gstrickd::Metadata *first_ledger(const std::vector<gstrickd::Record> &records) {
  for (const auto &record : records) {
    if (record.type == RecordType::ledger) {
      return &record.metadata;
    }
  }
  return nullptr;
}

std::string text_of(const gstrickd::Metadata &m, const std::string &key) {
  const auto it = m.find(key);
  if (it == m.end() || !std::holds_alternative<std::string>(it->second)) {
    return "<missing>";
  }
  return std::get<std::string>(it->second);
}

std::size_t array_size_of(const gstrickd::Metadata &m, const std::string &key) {
  const auto it = m.find(key);
  if (it == m.end() || !std::holds_alternative<std::vector<double>>(it->second)) {
    return 0;
  }
  return std::get<std::vector<double>>(it->second).size();
}

bool throws_invalid(const std::function<void()> &fn) {
  try {
    fn();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

// Ledger keys and record order: stream, run, ledger, event, ...
void test_ledger_layout() {
  const auto input = ground_stream(1000, 1000);
  auto config = config_for(zero_table());
  config.envelope = {constant_factor("qe", 0.9), constant_factor("mirror", 0.5)};
  const auto result = propagate(input, config);
  const auto records = records_of(result.bytes);
  CHECK(records.size() > 3);
  if (records.size() <= 3) {
    return;
  }
  CHECK(records[0].type == RecordType::stream);
  CHECK(records[1].type == RecordType::run);
  CHECK(records[2].type == RecordType::ledger);
  CHECK(records[3].type == RecordType::event);
  const auto *ledger = first_ledger(records);
  CHECK(ledger != nullptr);
  if (ledger == nullptr) {
    return;
  }
  CHECK(text_of(*ledger, "factors") == "atmosphere,wavelength,qe,mirror");
  CHECK(text_of(*ledger, "atmosphere.state") == "applied");
  CHECK(text_of(*ledger, "atmosphere.source") == "test table");
  CHECK(text_of(*ledger, "wavelength.state") == "applied");
  CHECK(text_of(*ledger, "qe.state") == "envelope");
  CHECK(text_of(*ledger, "qe.source") == "qe.dat");
  CHECK(array_size_of(*ledger, "qe.wavelength_nm") == 2);
  CHECK(array_size_of(*ledger, "mirror.value") == 2);
  CHECK(text_of(*ledger, "mirror.state") == "envelope");
  // Run keys: no profile, so refraction_applied is passed through unchanged.
  CHECK(std::get<std::int64_t>(records[1].metadata.at("refraction_applied")) == 0);
  CHECK(result.summary.photons_in == 2000);
  CHECK(result.summary.dropped == 0);
  // Without store_p the output has no p column.
  CHECK(photon_column(records, 1, Column::p).empty());
}

// Thinned output count against the weighted expectation (same config, same wavelengths).
void test_expected_count() {
  const auto input = ground_stream(60000, 20000);
  auto config = config_for(real_table());
  config.store_p = true;
  config.mode = dlufd::OutputMode::thinned;
  const auto thinned = propagate(input, config);
  config.mode = dlufd::OutputMode::weighted;
  const auto weighted = propagate(input, config);
  const auto thinned_records = records_of(thinned.bytes);
  const auto weighted_records = records_of(weighted.bytes);

  for (const std::uint32_t telescope : {1U, 2U}) {
    const auto w = photon_column(weighted_records, telescope, Column::w);
    const auto p = photon_column(weighted_records, telescope, Column::p);
    CHECK(w.size() == p.size());
    double mean = 0.0;     // sum w p
    double variance = 0.0; // sum w p (1 - p)
    for (std::size_t i = 0; i < w.size(); ++i) {
      mean += static_cast<double>(w[i]);
      variance += static_cast<double>(w[i]) * (1.0 - static_cast<double>(p[i]));
    }
    // w in the weighted output is w_in * p, so w_out (1 - p) = w_in p (1 - p).
    const double sigma = std::sqrt(variance);
    const auto count = static_cast<double>(photon_rows(thinned_records, telescope));
    CHECK(sigma > 100.0);
    CHECK_NEAR(count, mean, 4.0 * sigma);
    std::fprintf(stderr, "telescope %u: thinned %.0f, expected %.1f +- %.1f\n", telescope, count,
                 mean, sigma);
  }
  CHECK(weighted.summary.photons_out == 120000);
}

// Constant envelope 0.9 * 0.5 = 0.45 on a zero-absorption table: p = 0.45 exactly.
void test_weighted_analytic() {
  const auto input = ground_stream(5000, 5000);
  auto config = config_for(zero_table());
  config.envelope = {constant_factor("qe", 0.5), constant_factor("mirror", 0.8)};
  config.store_p = true;
  config.mode = dlufd::OutputMode::weighted;
  const auto weighted = propagate(input, config);
  // Sum of the input weights of both telescopes (ground_stream writes telescopes 1 and 2).
  double expected_weight = 0.0;
  for (std::size_t i = 0; i < 5000; ++i) {
    expected_weight += 2.0 * ground_row(i).w;
  }
  CHECK_NEAR(weighted.summary.weight_out, 0.4 * expected_weight, 1e-6 * expected_weight);
  const auto records = records_of(weighted.bytes);
  for (const float p : photon_column(records, 1, Column::p)) {
    CHECK_NEAR(static_cast<double>(p), 0.4, 1e-6);
  }

  // Thinned: Binomial with p = 0.4; variance sum w p (1 - p).
  config.mode = dlufd::OutputMode::thinned;
  const auto thinned = propagate(input, config);
  const double sigma = std::sqrt(0.24 * expected_weight);
  CHECK_NEAR(static_cast<double>(thinned.summary.photons_out), 0.4 * expected_weight, 4.0 * sigma);
}

// The input ledger must not show a factor applied already; pending or absent is fine.
void test_never_apply_twice() {
  const auto table = zero_table();
  const auto run_then_ledger = [](const gstrickd::Metadata &ledger) {
    return make_stream([&](gstrickd::Writer &writer) {
      writer.write_metadata(RecordType::run, run_keys(0));
      writer.write_metadata(RecordType::ledger, ledger);
      writer.write_metadata(RecordType::event, event_keys());
      writer.write_metadata(RecordType::event_end, {});
      writer.write_metadata(RecordType::run_end, {});
    });
  };
  gstrickd::Metadata applied;
  applied["factors"] = std::string("atmosphere");
  applied["atmosphere.state"] = std::string("applied");
  CHECK(throws_invalid([&] { propagate(run_then_ledger(applied), config_for(table)); }));

  gstrickd::Metadata qe_applied;
  qe_applied["factors"] = std::string("qe");
  qe_applied["qe.state"] = std::string("applied");
  auto config = config_for(table);
  config.envelope = {constant_factor("qe", 0.9)};
  CHECK(throws_invalid([&] { propagate(run_then_ledger(qe_applied), config); }));

  gstrickd::Metadata pending;
  pending["factors"] = std::string("atmosphere,qe");
  pending["atmosphere.state"] = std::string("pending");
  pending["qe.state"] = std::string("pending");
  bool accepted = true;
  try {
    propagate(run_then_ledger(pending), config);
  } catch (const std::exception &) {
    accepted = false;
  }
  CHECK(accepted);
}

// Output does not depend on how the input is split into chunks.
void test_chunk_size_determinism() {
  const auto table = real_table();
  auto config = config_for(table);
  config.store_p = true;
  const std::vector<std::string> ground = {ground_stream(50000, 1000), ground_stream(50000, 333),
                                           ground_stream(50000, gstrickd::max_chunk_rows)};
  std::vector<std::vector<float>> reference;
  for (std::size_t k = 0; k < ground.size(); ++k) {
    const auto records = records_of(propagate(ground[k], config).bytes);
    std::size_t index = 0;
    for (const std::uint32_t telescope : {1U, 2U}) {
      for (const Column c : {Column::x, Column::t, Column::w, Column::zem, Column::p}) {
        auto values = photon_column(records, telescope, c);
        if (k == 0) {
          reference.push_back(std::move(values));
        } else {
          CHECK(values == reference[index]);
        }
        ++index;
      }
    }
  }
  CHECK(!reference.empty() && !reference[0].empty());

  // Sky photons (some upwards) with chunk sizes 777 and 2^16.
  const auto sky_stream = [](std::size_t n, std::size_t chunk_rows) {
    return make_stream([&](gstrickd::Writer &writer) {
      writer.write_metadata(RecordType::run, run_keys(0));
      writer.write_metadata(RecordType::event, event_keys());
      writer.write_metadata(RecordType::telescope, telescope_keys(2, Frame::sky));
      for (std::size_t first = 0; first < n; first += chunk_rows) {
        const std::size_t rows = std::min(chunk_rows, n - first);
        writer.write_chunk(make_chunk(Frame::sky, 2, first, rows, sky_row));
      }
      writer.write_metadata(RecordType::telescope_end, telescope_end_keys(2));
      writer.write_metadata(RecordType::event_end, {});
      writer.write_metadata(RecordType::run_end, {});
    });
  };
  const auto small = propagate(sky_stream(20000, 777), config);
  const auto large = propagate(sky_stream(20000, gstrickd::max_chunk_rows), config);
  CHECK(small.summary.dropped == large.summary.dropped);
  CHECK(small.summary.dropped > 1000);
  const auto small_records = records_of(small.bytes);
  const auto large_records = records_of(large.bytes);
  for (const Column c : {Column::x, Column::y, Column::cx, Column::t, Column::zem, Column::w}) {
    CHECK(photon_column(small_records, 2, c) == photon_column(large_records, 2, c));
  }
  CHECK(!photon_column(small_records, 2, Column::x).empty());
}

// Analytic kinematics: 10 km above telescope 1 (at x = 100 m), direction (0.2, -0.1, -cz),
// telescope plane at z = 0, no profile. The photon arrives at the straight-line position
// and time t_in + s / c, with s the path length to the plane.
void test_sky_analytic() {
  const double cx = 0.2;
  const double cy = -0.1;
  const double cz_magnitude = std::sqrt(1.0 - cx * cx - cy * cy);
  const double s = 10000.0 / cz_magnitude;
  const auto input = make_stream([&](gstrickd::Writer &writer) {
    writer.write_metadata(RecordType::run, run_keys(-1));
    writer.write_metadata(RecordType::event, event_keys());
    writer.write_metadata(RecordType::telescope, telescope_keys(1, Frame::sky));
    const std::function<Row(std::size_t)> rows = [&](std::size_t i) {
      Row row;
      if (i == 0) {
        row.x = 150.0;
        row.y = -40.0;
        row.z = 10000.0;
        row.cx = cx;
        row.cy = cy;
        row.cz = -cz_magnitude;
        row.t = 5.0;
      } else if (i == 1) { // upwards
        row = Row{0.0, 0.0, 10000.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
      } else { // emitted below the telescope plane
        row = Row{0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0};
      }
      return row;
    };
    writer.write_chunk(make_chunk(Frame::sky, 1, 0, 3, rows));
    writer.write_metadata(RecordType::telescope_end, telescope_end_keys(1));
    writer.write_metadata(RecordType::event_end, {});
    writer.write_metadata(RecordType::run_end, {});
  });
  auto config = config_for(zero_table());
  const auto result = propagate(input, config);
  CHECK(result.summary.dropped == 2);
  const auto records = records_of(result.bytes);
  CHECK(photon_rows(records, 1) == 1);
  // Weighted output drops the dropped rows (w = 0) and keeps the arrived photon.
  config.mode = dlufd::OutputMode::weighted;
  CHECK(propagate(input, config).summary.photons_out == 1);
  const auto x = photon_column(records, 1, Column::x);
  const auto y = photon_column(records, 1, Column::y);
  const auto t = photon_column(records, 1, Column::t);
  const auto zem = photon_column(records, 1, Column::zem);
  const auto cz_column = photon_column(records, 1, Column::cz);
  CHECK(x.size() == 1 && y.size() == 1 && t.size() == 1 && zem.size() == 1);
  if (x.size() != 1 || y.size() != 1 || t.size() != 1 || zem.size() != 1) {
    return;
  }
  CHECK(cz_column.empty()); // frame G: z and cz are not written
  CHECK_NEAR(x[0], 150.0 + s * cx - 100.0, 1e-3);
  CHECK_NEAR(y[0], -40.0 + s * cy, 1e-3);
  CHECK_NEAR(t[0], 5.0 + s / dlufd::speed_of_light_m_per_ns, 1e-2);
  CHECK_NEAR(zem[0], 12000.0, 1e-3);
  // Frame of the output telescope record is G.
  bool ground_telescope = false;
  for (const auto &record : records) {
    if (record.type == RecordType::telescope) {
      ground_telescope = std::get<std::int64_t>(record.metadata.at("frame")) == 1;
    }
  }
  CHECK(ground_telescope);
}

// Refraction: the photon at 60 degrees zenith. Expected position and time are the refraction
// kernel's output (bend) at the telescope plane; without refraction the air delay is added.
void test_sky_refraction() {
  const double cx = std::sqrt(3.0) / 2.0; // cz = -1/2: 60 degrees zenith
  const double cz_magnitude = 0.5;
  const double s = 10000.0 / cz_magnitude;
  const auto profile = ecmwf_profile();
  const auto input = make_stream([&](gstrickd::Writer &writer) {
    writer.write_metadata(RecordType::run, run_keys(0));
    writer.write_metadata(RecordType::event, event_keys());
    writer.write_metadata(RecordType::telescope, telescope_keys(1, Frame::sky));
    writer.write_chunk(make_chunk(Frame::sky, 1, 0, 1, [&](std::size_t) {
      return Row{150.0, -40.0, 10000.0, cx, 0.0, -cz_magnitude, 5.0, 0.0, 1.0};
    }));
    writer.write_metadata(RecordType::telescope_end, telescope_end_keys(1));
    writer.write_metadata(RecordType::event_end, {});
    writer.write_metadata(RecordType::run_end, {});
  });
  auto config = config_for(zero_table());
  config.profile.emplace(profile);
  config.profile_source = "ecmwf";
  const auto result = propagate(input, config);
  const auto records = records_of(result.bytes);
  CHECK(photon_rows(records, 1) == 1);
  const auto x = photon_column(records, 1, Column::x);
  const auto y = photon_column(records, 1, Column::y);
  const auto t = photon_column(records, 1, Column::t);
  if (x.size() != 1 || t.size() != 1) {
    CHECK(false);
    return;
  }
  // Kernel reference: bend from the emission height to the plane at the site altitude.
  const dlufd::RefractionTables tables(profile, site_m);
  double dx = cx;
  double dy = 0.0;
  double xm = 150.0 + s * cx - 100.0;
  double ym = -40.0;
  double tb = 0.0;
  CHECK(dlufd::bend(tables, 12000.0, dx, dy, xm, ym, tb));
  const double analytic_x = 150.0 + s * cx - 100.0;
  CHECK(std::abs(xm - analytic_x) > 1.0); // about 3 m at 60 degrees
  // Columns are float32: at 17 km one ulp is about 2 mm.
  CHECK_NEAR(static_cast<double>(x[0]), xm, 5e-3);
  CHECK_NEAR(static_cast<double>(y[0]), ym, 1e-3);
  CHECK_NEAR(static_cast<double>(t[0]), 5.0 + tb, 2e-2);

  // The same photon with refraction already applied upstream: no bending, air delay added.
  auto applied_input = make_stream([&](gstrickd::Writer &writer) {
    writer.write_metadata(RecordType::run, run_keys(1));
    writer.write_metadata(RecordType::event, event_keys());
    writer.write_metadata(RecordType::telescope, telescope_keys(1, Frame::sky));
    writer.write_chunk(make_chunk(Frame::sky, 1, 0, 1, [&](std::size_t) {
      return Row{150.0, -40.0, 10000.0, cx, 0.0, -cz_magnitude, 5.0, 0.0, 1.0};
    }));
    writer.write_metadata(RecordType::telescope_end, telescope_end_keys(1));
    writer.write_metadata(RecordType::event_end, {});
    writer.write_metadata(RecordType::run_end, {});
  });
  const auto plain = records_of(propagate(applied_input, config).bytes);
  const double zenith = std::acos(cz_magnitude);
  const double air_ns = dlufd::straight_path_air_delay_ns(profile, site_m, 12000.0, zenith);
  CHECK_NEAR(static_cast<double>(photon_column(plain, 1, Column::x)[0]), analytic_x, 5e-3);
  CHECK_NEAR(static_cast<double>(photon_column(plain, 1, Column::t)[0]),
             5.0 + s / dlufd::speed_of_light_m_per_ns + air_ns, 2e-2);
  // The output run record now says refraction is applied.
  for (const auto &record : records) {
    if (record.type == RecordType::run) {
      CHECK(std::get<std::int64_t>(record.metadata.at("refraction_applied")) == 1);
    }
  }
}

// Sky photons need the z column; photons of an unknown telescope are rejected.
void test_rejects_invalid_input() {
  const auto no_z = make_stream([&](gstrickd::Writer &writer) {
    writer.write_metadata(RecordType::run, run_keys(0));
    writer.write_metadata(RecordType::event, event_keys());
    writer.write_metadata(RecordType::telescope, telescope_keys(1, Frame::sky));
    gstrickd::PhotonChunk chunk = make_chunk(Frame::sky, 1, 0, 2, sky_row);
    chunk.resize(2, gstrickd::required_columns);
    writer.write_chunk(chunk);
  });
  CHECK(throws_invalid([&] { propagate(no_z, config_for(zero_table())); }));

  const auto unknown = make_stream([&](gstrickd::Writer &writer) {
    writer.write_metadata(RecordType::run, run_keys(0));
    writer.write_metadata(RecordType::event, event_keys());
    writer.write_metadata(RecordType::telescope, telescope_keys(1, Frame::ground));
    writer.write_chunk(make_chunk(Frame::ground, 2, 0, 2, ground_row));
  });
  CHECK(throws_invalid([&] { propagate(unknown, config_for(zero_table())); }));
}

// ---- Single scattering -------------------------------------------------------------------

constexpr std::uint64_t scatter_seed = 20261107;

// Scattering table of the test data: 340 to 400 nm, exponential atmosphere (docs/scattering.md).
dlufd::ScatteringTable scattering_table() {
  return dlufd::ScatteringTable::from_file(std::string(DLUFD_TEST_DATA_DIR) +
                                           "/scattering_exponential.dat");
}

// Telescope 1 at 3 km east of the beam, telescope 2 at 2.5 km south; spheres of 12 m radius.
gstrickd::Metadata scatter_run_keys() {
  gstrickd::Metadata m;
  m["run_number"] = i64(1);
  m["observation_level_m"] = site_m;
  m["telescope_x_m"] = std::vector<double>{3000.0, 0.0};
  m["telescope_y_m"] = std::vector<double>{0.0, -2500.0};
  m["telescope_z_m"] = std::vector<double>{0.0, 0.0};
  m["telescope_r_m"] = std::vector<double>{12.0, 12.0};
  m["refraction_applied"] = i64(0);
  m["wavelength_min_nm"] = 250.0;
  m["wavelength_max_nm"] = 700.0;
  return m;
}

// Vertical beam (cx = cy = 0) at x = offset_x, emitted 15 km above the observation level,
// wavelength 355 nm, weights 1 to 3.
gstrickd::PhotonChunk beam_chunk(std::uint32_t telescope, std::size_t first, std::size_t n,
                                 double offset_x) {
  gstrickd::PhotonChunk chunk;
  chunk.telescope_id = telescope;
  chunk.frame = Frame::sky;
  chunk.first_index = first;
  const std::uint32_t mask = gstrickd::column_bit(Column::x) | gstrickd::column_bit(Column::y) |
                             gstrickd::column_bit(Column::z) | gstrickd::column_bit(Column::cx) |
                             gstrickd::column_bit(Column::cy) | gstrickd::column_bit(Column::cz) |
                             gstrickd::column_bit(Column::t) | gstrickd::column_bit(Column::w) |
                             gstrickd::column_bit(Column::wavelength);
  chunk.resize(n, mask);
  for (std::size_t r = 0; r < n; ++r) {
    const std::size_t i = first + r;
    chunk[Column::x][r] = static_cast<float>(offset_x);
    chunk[Column::y][r] = 0.0F;
    chunk[Column::z][r] = 15000.0F;
    chunk[Column::cx][r] = 0.0F;
    chunk[Column::cy][r] = 0.0F;
    chunk[Column::cz][r] = -1.0F;
    chunk[Column::t][r] = 0.0F;
    chunk[Column::w][r] = static_cast<float>(1 + i % 3);
    chunk[Column::wavelength][r] = 355.0F;
  }
  return chunk;
}

// One event. With two_blocks, telescope 2 has its own beam (at x = 200 m) and primaries; else
// telescope 2 has no block and receives its scattered photons in a block of its own.
std::string beam_stream(std::size_t n, std::size_t chunk_rows, bool two_blocks,
                        gstrickd::Metadata run = scatter_run_keys()) {
  return make_stream([&](gstrickd::Writer &writer) {
    writer.write_metadata(RecordType::run, run);
    writer.write_metadata(RecordType::event, event_keys());
    writer.write_metadata(RecordType::telescope, telescope_keys(1, Frame::sky));
    for (std::size_t first = 0; first < n; first += chunk_rows) {
      writer.write_chunk(beam_chunk(1, first, std::min(chunk_rows, n - first), 0.0));
    }
    writer.write_metadata(RecordType::telescope_end, telescope_end_keys(1));
    if (two_blocks) {
      writer.write_metadata(RecordType::telescope, telescope_keys(2, Frame::sky));
      for (std::size_t first = 0; first < n; first += chunk_rows) {
        writer.write_chunk(beam_chunk(2, first, std::min(chunk_rows, n - first), 200.0));
      }
      writer.write_metadata(RecordType::telescope_end, telescope_end_keys(2));
    }
    writer.write_metadata(RecordType::event_end, {});
    writer.write_metadata(RecordType::run_end, {});
  });
}

dlufd::PropagatorConfig scattering_config(const dlufd::TransmissionTable &table, bool weighted) {
  auto config = config_for(table);
  config.mode = weighted ? dlufd::OutputMode::weighted : dlufd::OutputMode::thinned;
  config.envelope.push_back(constant_factor("qe", 0.5));
  config.scattering.enabled = true;
  config.scattering.points_per_photon = 8;
  config.scattering_table.emplace(scattering_table());
  config.scattering_source = "scattering_exponential.dat";
  return config;
}

// Sum of the telescope_end photons (sum of w) of one telescope; 0 without a block.
double telescope_photons(const std::vector<gstrickd::Record> &records, std::uint32_t telescope) {
  double total = 0.0;
  for (const auto &record : records) {
    if (record.type == RecordType::telescope_end &&
        std::get<std::int64_t>(record.metadata.at("telescope_id")) == telescope) {
      total += std::get<double>(record.metadata.at("photons"));
    }
  }
  return total;
}

// first_index of the photon chunks of one telescope counts its rows from 0 without gaps.
bool contiguous_indices(const std::vector<gstrickd::Record> &records, std::uint32_t telescope) {
  std::uint64_t next = 0;
  for (const auto &record : records) {
    if (record.type == RecordType::photons && record.chunk.telescope_id == telescope) {
      if (record.chunk.first_index != next) {
        return false;
      }
      next += record.chunk.rows();
    }
  }
  return true;
}

// Scattered weight per telescope from scatter_chunk called directly on the n beam photons.
std::vector<double> kernel_weights(const dlufd::TransmissionTable &table, std::size_t n) {
  dlufd::ScatteringConfig config;
  config.enabled = true;
  config.points_per_photon = 8;
  const std::vector<dlufd::ScatterTelescope> sphere = {{1, 3000.0, 0.0, 0.0, 12.0},
                                                       {2, 0.0, -2500.0, 0.0, 12.0}};
  std::vector<std::vector<gstrickd::PhotonChunk>> lists(sphere.size());
  dlufd::scatter_chunk(config, scattering_table(), table, sphere, beam_chunk(1, 0, n, 0.0),
                       scatter_seed, 7, 0, lists);
  std::vector<double> sums(sphere.size(), 0.0);
  for (std::size_t i = 0; i < sphere.size(); ++i) {
    for (const auto &chunk : lists[i]) {
      for (const float w : chunk[Column::w]) {
        sums[i] += static_cast<double>(w);
      }
    }
  }
  return sums;
}

// Scattered photons equal the kernel result times the envelope (E = 0.5), added to the block of
// their telescope; the primaries are unchanged; first_index stays contiguous.
void test_scattering_weights() {
  const auto table = real_table();
  const std::size_t n = 400;
  const std::string input = beam_stream(n, n, false);
  auto off = scattering_config(table, true);
  off.scattering = dlufd::ScatteringConfig{};
  off.scattering_table.reset();
  const auto without = records_of(propagate(input, off).bytes);
  const auto with = records_of(propagate(input, scattering_config(table, true)).bytes);
  const auto expected = kernel_weights(table, n);
  CHECK(expected[0] > 0.0 && expected[1] > 0.0);
  for (const std::uint32_t telescope : {1U, 2U}) {
    const double added = telescope_photons(with, telescope) - telescope_photons(without, telescope);
    CHECK_NEAR(added, 0.5 * expected[telescope - 1], 1e-5 * expected[telescope - 1]);
    CHECK(contiguous_indices(with, telescope));
  }
  CHECK(photon_rows(with, 2) > 0); // block written for telescope 2 without input block
  CHECK(photon_rows(without, 2) == 0);
  CHECK(photon_column(with, 1, Column::x).size() > photon_column(without, 1, Column::x).size());
}

// Extinction disabled by an all-zero transmission table: the scattered weights are those of
// the kernel with that table times the envelope. The path transmission T(P -> C) of the kernel
// uses the same table (it is in w), so it is not compared with the real-table run.
void test_scattering_not_attenuated_twice() {
  const auto zero = zero_table();
  const std::size_t n = 300;
  const std::string input = beam_stream(n, n, false); // one block: the kernel reference is one beam
  auto off = scattering_config(zero, true);
  off.scattering = dlufd::ScatteringConfig{};
  off.scattering_table.reset();
  const auto without = records_of(propagate(input, off).bytes);
  const auto with = records_of(propagate(input, scattering_config(zero, true)).bytes);
  const auto expected = kernel_weights(zero, n);
  for (const std::uint32_t telescope : {1U, 2U}) {
    const double added = telescope_photons(with, telescope) - telescope_photons(without, telescope);
    CHECK_NEAR(added, 0.5 * expected[telescope - 1], 1e-5 * expected[telescope - 1]);
  }
}

void test_scattering_ledger() {
  const auto table = real_table();
  const std::string input = beam_stream(100, 100, false);
  const auto on = records_of(propagate(input, scattering_config(table, false)).bytes);
  const auto *ledger = first_ledger(on);
  CHECK(ledger != nullptr);
  CHECK(text_of(*ledger, "factors") == "atmosphere,wavelength,scattering,qe");
  CHECK(text_of(*ledger, "scattering.state") == "applied");
  CHECK(text_of(*ledger, "scattering.source") == "scattering_exponential.dat");
  CHECK(std::get<std::int64_t>(ledger->at("scattering.points_per_photon")) == 8);

  auto off = scattering_config(table, false);
  off.scattering = dlufd::ScatteringConfig{};
  off.scattering_table.reset();
  const auto plain_records = records_of(propagate(input, off).bytes);
  const auto *plain = first_ledger(plain_records);
  CHECK(plain != nullptr);
  CHECK(plain == nullptr || text_of(*plain, "factors") == "atmosphere,wavelength,qe");
  CHECK(plain == nullptr || plain->count("scattering.state") == 0);
  CHECK(plain == nullptr || plain->count("scattering.points_per_photon") == 0);
}

// Thinned and weighted output do not depend on how the input is split into chunks, with
// primaries in two telescope blocks (the scatter index runs over both).
void test_scattering_chunk_determinism() {
  const auto table = real_table();
  for (const bool weighted : {false, true}) {
    const auto config = scattering_config(table, weighted);
    const auto reference = records_of(propagate(beam_stream(3000, 3000, true), config).bytes);
    CHECK(photon_rows(reference, 2) > 0);
    for (const std::size_t chunk_rows : {std::size_t{1000}, std::size_t{333}}) {
      const auto records = records_of(propagate(beam_stream(3000, chunk_rows, true), config).bytes);
      for (const std::uint32_t telescope : {1U, 2U}) {
        for (const Column c : {Column::x, Column::y, Column::t, Column::w, Column::zem}) {
          CHECK(photon_column(records, telescope, c) == photon_column(reference, telescope, c));
        }
        CHECK(contiguous_indices(records, telescope));
      }
    }
  }
}

void test_scattering_rejects() {
  const auto table = real_table();
  const std::string input = beam_stream(10, 10, false);
  auto no_table = scattering_config(table, false);
  no_table.scattering_table.reset();
  CHECK(throws_invalid([&] { propagate(input, no_table); }));

  gstrickd::Metadata no_radius = scatter_run_keys();
  no_radius.erase("telescope_r_m");
  CHECK(throws_invalid(
      [&] { propagate(beam_stream(10, 10, false, no_radius), scattering_config(table, false)); }));

  const std::string ground = ground_stream(100, 100);
  CHECK(throws_invalid([&] { propagate(ground, scattering_config(table, false)); }));
}

} // namespace

int main() {
  return run_checks([] {
    test_ledger_layout();
    test_expected_count();
    test_weighted_analytic();
    test_never_apply_twice();
    test_chunk_size_determinism();
    test_sky_analytic();
    test_sky_refraction();
    test_rejects_invalid_input();
    test_scattering_weights();
    test_scattering_not_attenuated_twice();
    test_scattering_ledger();
    test_scattering_chunk_determinism();
    test_scattering_rejects();
  });
}
