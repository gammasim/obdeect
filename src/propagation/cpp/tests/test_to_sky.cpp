#include "check.hpp"

#include "obdeect/propagation/ground_to_sky.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/propagator.hpp"
#include "obdeect/propagation/timing.hpp"

#include "obdeect/stream/constants.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

using obdeect::propagation::AtmosphereProfile;
using obdeect::propagation::ground_slice;
using obdeect::propagation::ground_to_sky;
using obdeect::propagation::ground_to_sky_row;
using obdeect::propagation::ground_to_sky_stream;
using obdeect::propagation::GroundPlane;
using obdeect::stream::Column;
using obdeect::stream::Frame;
using obdeect::stream::PhotonChunk;
using obdeect::stream::RecordType;

constexpr double site_m = 2000.0;
constexpr double c_m_per_ns = obdeect::stream::speed_of_light_m_per_ns;

AtmosphereProfile ecmwf_profile() {
  return AtmosphereProfile::from_file(std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) +
                                      "/atmprof_ecmwf_north_winter_fixed.dat");
}

// Two telescopes: (100 m, 0 m, 0 m) and (-50 m, 300 m, 25 m) relative to the array centre.
obdeect::stream::Metadata run_keys(int refraction) {
  obdeect::stream::Metadata m;
  m["run_number"] = std::int64_t{1};
  m["observation_level_m"] = site_m;
  m["telescope_x_m"] = std::vector<double>{100.0, -50.0};
  m["telescope_y_m"] = std::vector<double>{0.0, 300.0};
  m["telescope_z_m"] = std::vector<double>{0.0, 25.0};
  m["telescope_r_m"] = std::vector<double>{12.0, 12.0};
  if (refraction >= 0) {
    m["refraction_applied"] = std::int64_t{refraction};
  }
  m["wavelength_min_nm"] = 250.0;
  m["wavelength_max_nm"] = 700.0;
  return m;
}

// One photon. Sky photons use z; G photons use zem. cz < -2 means "no cz column".
struct Row {
  double x = 0.0, y = 0.0, z = 0.0, cx = 0.0, cy = 0.0, cz = -2.0, t = 0.0, zem = 0.0, w = 1.0;
};

// Chunk of frame G (x y cx cy t w zem, cz optional) or sky (x y z cx cy t w, cz optional).
PhotonChunk make_chunk(Frame frame, std::uint32_t telescope, std::uint64_t first, std::size_t n,
                       const std::function<Row(std::size_t)> &row_of, bool with_cz) {
  PhotonChunk chunk;
  chunk.telescope_id = telescope;
  chunk.frame = frame;
  chunk.first_index = first;
  const bool sky = frame == Frame::sky;
  std::uint32_t mask =
      obdeect::stream::column_bit(Column::x) | obdeect::stream::column_bit(Column::y) |
      obdeect::stream::column_bit(Column::cx) | obdeect::stream::column_bit(Column::cy) |
      obdeect::stream::column_bit(Column::t) | obdeect::stream::column_bit(Column::w);
  mask |= sky ? obdeect::stream::column_bit(Column::z) : obdeect::stream::column_bit(Column::zem);
  if (with_cz) {
    mask |= obdeect::stream::column_bit(Column::cz);
  }
  chunk.resize(n, mask);
  for (std::size_t r = 0; r < n; ++r) {
    const Row row = row_of(r);
    chunk[Column::x][r] = static_cast<float>(row.x);
    chunk[Column::y][r] = static_cast<float>(row.y);
    chunk[Column::cx][r] = static_cast<float>(row.cx);
    chunk[Column::cy][r] = static_cast<float>(row.cy);
    chunk[Column::t][r] = static_cast<float>(row.t);
    chunk[Column::w][r] = static_cast<float>(row.w);
    if (sky) {
      chunk[Column::z][r] = static_cast<float>(row.z);
    } else {
      chunk[Column::zem][r] = static_cast<float>(row.zem);
    }
    if (with_cz) {
      chunk[Column::cz][r] = static_cast<float>(row.cz);
    }
  }
  return chunk;
}

// Sky photons 3 to 23 km above the observation level, all moving downwards.
Row sky_row(std::size_t i) {
  Row row;
  row.x = 500.0 * std::sin(0.1 * static_cast<double>(i));
  row.y = 500.0 * std::cos(0.07 * static_cast<double>(i));
  row.z = 3000.0 + 5000.0 * static_cast<double>(i % 5);
  row.cx = 0.3 * std::sin(0.01 * static_cast<double>(i));
  row.cy = 0.3 * std::cos(0.013 * static_cast<double>(i));
  row.cz = -std::sqrt(1.0 - row.cx * row.cx - row.cy * row.cy);
  row.t = static_cast<double>(i % 11); // small sky times keep float32 rounding of G times exact
  row.w = 1.0 + static_cast<double>(i % 3);
  return row;
}

// G photons on a spread of positions, directions and emission altitudes (cz from the direction).
Row ground_row(std::size_t i) {
  Row row;
  row.x = 300.0 * std::sin(0.37 * static_cast<double>(i));
  row.y = 300.0 * std::cos(0.23 * static_cast<double>(i));
  row.cx = 0.4 * std::sin(0.001 * static_cast<double>(i));
  row.cy = 0.4 * std::cos(0.0013 * static_cast<double>(i));
  row.cz = -std::sqrt(1.0 - row.cx * row.cx - row.cy * row.cy);
  row.t = 30000.0 + static_cast<double>(i % 13);
  row.zem = 2100.0 + 1000.0 * static_cast<double>(i % 7); // above the plane at 2000 m
  row.w = 1.0 + static_cast<double>(i % 4);
  return row;
}

struct TelescopeData {
  std::uint32_t id;
  Frame frame;
  obdeect::stream::Metadata keys; // extra telescope record keys (x_m, ...)
  std::vector<PhotonChunk> chunks;
};

std::string build_stream(const obdeect::stream::Metadata &run,
                         const std::vector<TelescopeData> &telescopes) {
  std::ostringstream out;
  {
    obdeect::stream::Writer writer(out, "test");
    writer.write_metadata(RecordType::run, run);
    obdeect::stream::Metadata event;
    event["event_number"] = std::int64_t{7};
    writer.write_metadata(RecordType::event, event);
    for (const auto &telescope : telescopes) {
      obdeect::stream::Metadata keys = telescope.keys;
      keys["telescope_id"] = static_cast<std::int64_t>(telescope.id);
      keys["frame"] = static_cast<std::int64_t>(telescope.frame);
      writer.write_metadata(RecordType::telescope, keys);
      for (PhotonChunk chunk : telescope.chunks) {
        chunk.telescope_id = telescope.id;
        writer.write_chunk(chunk);
      }
      obdeect::stream::Metadata end;
      end["telescope_id"] = static_cast<std::int64_t>(telescope.id);
      end["photons"] = 0.0;
      writer.write_metadata(RecordType::telescope_end, end);
    }
    writer.write_metadata(RecordType::event_end, obdeect::stream::Metadata{});
    writer.write_metadata(RecordType::run_end, obdeect::stream::Metadata{});
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

std::size_t photon_rows(const std::vector<obdeect::stream::Record> &records,
                        std::uint32_t telescope) {
  std::size_t rows = 0;
  for (const auto &record : records) {
    if (record.type == RecordType::photons && record.chunk.telescope_id == telescope) {
      rows += record.chunk.rows();
    }
  }
  return rows;
}

// True if f throws std::invalid_argument.
template <class F> bool rejects(F f) {
  try {
    f();
  } catch (const std::invalid_argument &) {
    return true;
  } catch (...) {
    return false;
  }
  return false;
}

// Output of the stage (geometry only, no extinction) for a sky stream.
std::string stage_geometry(const std::string &input,
                           const std::optional<AtmosphereProfile> &profile) {
  obdeect::propagation::PropagatorConfig config;
  config.geometry_only = true;
  config.mode = obdeect::propagation::OutputMode::weighted;
  config.site_altitude_m = site_m;
  config.seed = 20261120;
  config.profile = profile;
  std::istringstream in(input);
  std::ostringstream out;
  obdeect::propagation::propagate_stream(in, out, config);
  return out.str();
}

// Analytic cases: the row function in double precision.
void analytic_cases() {
  GroundPlane plane{100.0, -50.0, 0.0, site_m, nullptr};
  double x = 0.0, y = 0.0, z = 0.0, t = 0.0;

  // Vertical photon, 300 m above the plane: s = 300 m, sky position directly above.
  CHECK(ground_to_sky_row(plane, 3.0, 4.0, 0.0, 0.0, site_m + 300.0, 5.0, x, y, z, t));
  CHECK_NEAR(x, 103.0, 1e-9);
  CHECK_NEAR(y, -46.0, 1e-9);
  CHECK_NEAR(z, 300.0, 1e-9);
  CHECK_NEAR(t, 5.0 - 300.0 / c_m_per_ns, 1e-9);

  // Inclined photon: cx = 0.6, cz = -0.8, 800 m above the plane, so s = 1000 m.
  CHECK(ground_to_sky_row(plane, 3.0, 4.0, 0.6, 0.0, site_m + 800.0, 10.0, x, y, z, t));
  CHECK_NEAR(x, 100.0 + 3.0 - 1000.0 * 0.6, 1e-9);
  CHECK_NEAR(y, -50.0 + 4.0, 1e-9);
  CHECK_NEAR(z, 800.0, 1e-9);
  CHECK_NEAR(t, 10.0 - 1000.0 / c_m_per_ns, 1e-9);

  // Telescope offset in z: h_t = 2025 m, so the path starts at the telescope plane.
  GroundPlane raised{0.0, 0.0, 25.0, site_m, nullptr};
  CHECK(ground_to_sky_row(raised, 1.0, 2.0, 0.0, 0.0, 2100.0, 0.0, x, y, z, t));
  CHECK_NEAR(x, 1.0, 1e-9);
  CHECK_NEAR(y, 2.0, 1e-9);
  CHECK_NEAR(z, 100.0, 1e-9);
  CHECK_NEAR(t, -75.0 / c_m_per_ns, 1e-9);
  // zem at the plane (zem <= h_t) and above the plane but horizontal: not convertible.
  CHECK(!ground_to_sky_row(raised, 1.0, 2.0, 0.0, 0.0, 2025.0, 0.0, x, y, z, t));
  CHECK(!ground_to_sky_row(plane, 1.0, 2.0, 1.0, 0.0, site_m + 500.0, 0.0, x, y, z, t));
  CHECK(!ground_to_sky_row(plane, std::nan(""), 2.0, 0.0, 0.0, site_m + 5.0, 0.0, x, y, z, t));
  CHECK(rejects(
      [&] { ground_to_sky_row(plane, 1.0, 2.0, 0.9, 0.9, site_m + 5.0, 0.0, x, y, z, t); }));
}

// Time inversion with a profile: the stage's straight-path air delay is subtracted.
void time_with_profile() {
  const AtmosphereProfile profile = ecmwf_profile();
  GroundPlane plane{100.0, -50.0, 0.0, site_m, &profile};
  double x = 0.0, y = 0.0, z = 0.0, t = 0.0;
  const double cx = 0.6;
  const double cz_magnitude = 0.8;
  const double zem = site_m + 800.0;
  CHECK(ground_to_sky_row(plane, 3.0, 4.0, cx, 0.0, zem, 10.0, x, y, z, t));
  const double zenith = std::atan2(cx, cz_magnitude);
  const double delay_ns = 1000.0 / c_m_per_ns + obdeect::propagation::straight_path_air_delay_ns(
                                                    profile, site_m, zem, zenith);
  CHECK(delay_ns > 1000.0 / c_m_per_ns);
  CHECK_NEAR(t, 10.0 - delay_ns, 1e-9);
  CHECK_NEAR(x, 100.0 + 3.0 - 1000.0 * cx, 1e-9); // the air delay does not move the position
}

// Chunk level: dropped rows are removed, counted, and the indices of the kept rows are kept.
void dropped_rows() {
  const std::uint64_t first = 100;
  const auto row_of = [](std::size_t r) {
    Row row = ground_row(r);
    row.cz = -std::sqrt(1.0 - row.cx * row.cx - row.cy * row.cy);
    if (r == 1) {
      row.zem = site_m; // at the observation level: below the plane
    }
    if (r == 3) {
      row.x = std::nan("");
    }
    if (r == 4) {
      row.cz = 0.5; // upwards, cz column present
    }
    if (r == 6) {
      row.cx = 1.0; // horizontal
      row.cy = 0.0;
      row.cz = 0.0;
    }
    return row;
  };
  const PhotonChunk input = make_chunk(Frame::ground, 2, first, 7, row_of, true);
  std::uint64_t dropped = 0;
  const auto pieces =
      ground_to_sky(input, GroundPlane{-50.0, 300.0, 0.0, site_m, nullptr}, dropped);
  CHECK(dropped == 4);
  CHECK(pieces.size() == 3);
  if (pieces.size() == 3) {
    CHECK(pieces[0].first_index == first && pieces[0].rows() == 1);
    CHECK(pieces[1].first_index == first + 2 && pieces[1].rows() == 1);
    CHECK(pieces[2].first_index == first + 5 && pieces[2].rows() == 1);
  }
  for (const auto &piece : pieces) {
    CHECK(piece.frame == Frame::sky);
    CHECK(piece.telescope_id == 2);
    CHECK(piece.has(Column::z));
    CHECK(piece.has(Column::zem));
    CHECK(piece.has(Column::cz));
  }
}

// Results depend only on the photon index: any chunking gives the same rows.
void chunk_size_independence() {
  const auto row_of = [](std::size_t r) {
    Row row = ground_row(r);
    if (r % 7 == 3) {
      row.zem = 1900.0; // below the plane: dropped
    }
    return row;
  };
  const std::size_t n = 300;
  const PhotonChunk whole =
      make_chunk(Frame::ground, 1, 1000, n, [&](std::size_t r) { return row_of(r); }, true);
  const GroundPlane plane{100.0, 0.0, 0.0, site_m, nullptr};

  using Output = std::array<float, 8>;
  const auto collect = [](const std::vector<PhotonChunk> &chunks) {
    std::map<std::uint64_t, Output> rows;
    for (const auto &chunk : chunks) {
      for (std::size_t r = 0; r < chunk.rows(); ++r) {
        const Output out{chunk[Column::x][r],   chunk[Column::y][r],  chunk[Column::z][r],
                         chunk[Column::cx][r],  chunk[Column::cy][r], chunk[Column::t][r],
                         chunk[Column::zem][r], chunk[Column::w][r]};
        CHECK(rows.count(chunk.first_index + r) == 0);
        rows[chunk.first_index + r] = out;
      }
    }
    return rows;
  };

  std::uint64_t dropped_whole = 0;
  const auto reference = collect(ground_to_sky(whole, plane, dropped_whole));
  CHECK(dropped_whole == n / 7 + 1); // rows 3, 10, ..., 297: 43 rows
  for (const std::size_t size : {std::size_t{1}, std::size_t{5}, std::size_t{64}}) {
    std::vector<PhotonChunk> pieces;
    std::uint64_t dropped = 0;
    for (std::size_t begin = 0; begin < n; begin += size) {
      const std::size_t end = std::min(n, begin + size);
      for (auto &piece : ground_to_sky(ground_slice(whole, begin, end), plane, dropped)) {
        pieces.push_back(piece);
      }
    }
    CHECK(dropped == dropped_whole);
    const auto rows = collect(pieces);
    CHECK(rows.size() == reference.size());
    for (const auto &[index, values] : reference) {
      const auto found = rows.find(index);
      CHECK(found != rows.end());
      if (found != rows.end()) {
        CHECK(found->second == values);
      }
    }
  }

  // Row order: the row function gives the same result for each photon in reverse order.
  std::vector<std::array<double, 4>> forward;
  for (std::size_t r = 0; r < n; ++r) {
    const Row row = row_of(r);
    std::array<double, 4> out{};
    if (!ground_to_sky_row(plane, row.x, row.y, row.cx, row.cy, row.zem, row.t, out[0], out[1],
                           out[2], out[3])) {
      out = {std::nan(""), 0.0, 0.0, 0.0};
    }
    forward.push_back(out);
  }
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t r = n - 1 - k;
    const Row row = row_of(r);
    std::array<double, 4> out{};
    if (!ground_to_sky_row(plane, row.x, row.y, row.cx, row.cy, row.zem, row.t, out[0], out[1],
                           out[2], out[3])) {
      out = {std::nan(""), 0.0, 0.0, 0.0};
    }
    CHECK(std::isnan(forward[r][0]) == std::isnan(out[0]));
    if (!std::isnan(out[0])) {
      CHECK(forward[r] == out);
    }
  }
}

// Invalid input is rejected with std::invalid_argument; the run and telescope records decide.
void rejected_inputs() {
  const PhotonChunk good = make_chunk(Frame::ground, 1, 0, 3, ground_row, true);
  const auto ok_stream = [&](const obdeect::stream::Metadata &run) {
    return build_stream(run, {TelescopeData{1, Frame::ground, {}, {good}}});
  };
  const auto convert = [](const std::string &bytes, const std::optional<AtmosphereProfile> &p) {
    std::istringstream in(bytes);
    std::ostringstream out;
    return ground_to_sky_stream(in, out, p);
  };

  // Curved input.
  obdeect::stream::Metadata curved = run_keys(0);
  curved["curvature_applied"] = std::int64_t{1};
  CHECK(rejects([&] { convert(ok_stream(curved), std::nullopt); }));

  // Missing observation level.
  obdeect::stream::Metadata no_level = run_keys(0);
  no_level.erase("observation_level_m");
  CHECK(rejects([&] { convert(ok_stream(no_level), std::nullopt); }));

  // Missing telescope positions.
  obdeect::stream::Metadata no_positions = run_keys(0);
  no_positions.erase("telescope_y_m");
  CHECK(rejects([&] { convert(ok_stream(no_positions), std::nullopt); }));

  // Telescope record in frame sky (or T): not G.
  CHECK(rejects([&] {
    convert(build_stream(run_keys(0), {TelescopeData{1, Frame::sky, {}, {good}}}), std::nullopt);
  }));

  // Telescope record position that differs from the run's position.
  obdeect::stream::Metadata wrong_position;
  wrong_position["x_m"] = 999.0;
  CHECK(rejects([&] {
    convert(build_stream(run_keys(0), {TelescopeData{1, Frame::ground, wrong_position, {good}}}),
            std::nullopt);
  }));
  // Telescope record whose position agrees is accepted.
  obdeect::stream::Metadata right_position;
  right_position["x_m"] = 100.0;
  CHECK(!rejects([&] {
    convert(build_stream(run_keys(0), {TelescopeData{1, Frame::ground, right_position, {good}}}),
            std::nullopt);
  }));

  // Telescope id outside the run.
  CHECK(rejects([&] {
    convert(build_stream(run_keys(0), {TelescopeData{3, Frame::ground, {}, {good}}}), std::nullopt);
  }));

  // A profile needs refraction_applied = 1 in the run.
  const auto profile = ecmwf_profile();
  CHECK(rejects([&] { convert(ok_stream(run_keys(0)), profile); }));
  CHECK(rejects([&] { convert(ok_stream(run_keys(-1)), profile); }));
  CHECK(!rejects([&] { convert(ok_stream(run_keys(1)), profile); }));
  CHECK(!rejects([&] { convert(ok_stream(run_keys(0)), std::nullopt); }));

  // Sky-frame photons and photons without zem are rejected at chunk level.
  std::uint64_t dropped = 0;
  const GroundPlane plane{100.0, -50.0, 0.0, site_m, nullptr};
  const PhotonChunk sky = make_chunk(Frame::sky, 1, 0, 2, sky_row, true);
  CHECK(rejects([&] { ground_to_sky(sky, plane, dropped); }));
  PhotonChunk no_zem = make_chunk(Frame::ground, 1, 0, 2, ground_row, false);
  no_zem[Column::zem].clear();
  CHECK(rejects([&] { ground_to_sky(no_zem, plane, dropped); }));
  // A cz column that does not match the direction is rejected.
  Row inconsistent = ground_row(0);
  inconsistent.cz = -0.1;
  const PhotonChunk bad_cz =
      make_chunk(Frame::ground, 1, 0, 1, [&](std::size_t) { return inconsistent; }, true);
  CHECK(rejects([&] { ground_to_sky(bad_cz, plane, dropped); }));
  CHECK(dropped == 0);
}

// Round trip through the stage: sky -> G (stage) -> sky (to-sky) -> G (stage), with and without a
// profile (refraction_applied = 1 with the profile, as CORSIKA writes it with refraction).
void round_trip(bool with_profile) {
  std::optional<AtmosphereProfile> profile;
  if (with_profile) {
    profile = ecmwf_profile();
  }
  const int refraction = with_profile ? 1 : 0;
  std::vector<TelescopeData> sky_telescopes;
  for (std::uint32_t id : {std::uint32_t{1}, std::uint32_t{2}}) {
    sky_telescopes.push_back(
        TelescopeData{id, Frame::sky, {}, {make_chunk(Frame::sky, id, 0, 400, sky_row, true)}});
  }
  const std::string sky_stream = build_stream(run_keys(refraction), sky_telescopes);

  // Sky -> G with the stage, then G -> sky with to-sky, then G again with the stage.
  const std::string ground = stage_geometry(sky_stream, profile);
  std::istringstream ground_in(ground);
  std::ostringstream back_out;
  const auto summary = ground_to_sky_stream(ground_in, back_out, profile);
  CHECK(summary.dropped == 0);
  CHECK(summary.photons_in == summary.photons_out);
  const std::string back_to_sky = back_out.str();
  const std::string ground_again = stage_geometry(back_to_sky, profile);

  const auto g1 = records_of(ground);
  const auto g2 = records_of(ground_again);
  const auto s2 = records_of(back_to_sky);
  for (const auto &record : s2) {
    if (record.type == RecordType::telescope) {
      CHECK(std::get<std::int64_t>(record.metadata.at("frame")) ==
            static_cast<std::int64_t>(Frame::sky));
    }
  }
  for (const std::uint32_t id : {std::uint32_t{1}, std::uint32_t{2}}) {
    CHECK(photon_rows(g1, id) == 400);
    CHECK(photon_rows(g2, id) == photon_rows(g1, id));
    for (const Column c : {Column::x, Column::y, Column::cx, Column::cy, Column::t, Column::zem}) {
      const auto a = photon_column(g1, id, c);
      const auto b = photon_column(g2, id, c);
      CHECK(a.size() == b.size());
      if (a.size() != b.size()) {
        continue;
      }
      const double tolerance = (c == Column::cx || c == Column::cy) ? 1e-6 : 1e-3;
      double worst = 0.0;
      for (std::size_t r = 0; r < a.size(); ++r) {
        worst = std::max(worst, std::abs(static_cast<double>(a[r]) - b[r]));
      }
      CHECK(worst <= tolerance);
    }
    CHECK(photon_column(g1, id, Column::w) == photon_column(g2, id, Column::w));
    // The back-converted sky photons are the stage's sky input (x, y, z to 1e-3 m, cx, cy to 1e-6).
    const PhotonChunk &original = sky_telescopes[id - 1].chunks[0];
    for (const Column c : {Column::x, Column::y, Column::z, Column::cx, Column::cy}) {
      const auto a = photon_column(s2, id, c);
      CHECK(a.size() == original.rows());
      if (a.size() != original.rows()) {
        continue;
      }
      const double tolerance = (c == Column::cx || c == Column::cy) ? 1e-6 : 1e-3;
      double worst = 0.0;
      for (std::size_t r = 0; r < a.size(); ++r) {
        worst = std::max(worst, std::abs(static_cast<double>(a[r]) - original[c][r]));
      }
      CHECK(worst <= tolerance);
    }
  }
}

} // namespace

int main() {
  return run_checks([] {
    analytic_cases();
    time_with_profile();
    dropped_rows();
    chunk_size_independence();
    rejected_inputs();
    round_trip(false);
    round_trip(true);
  });
}
