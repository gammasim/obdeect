#include "check.hpp"

#include "gstrickd/grisu.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace gstrickd;

namespace {

template <class Action> bool throws_error(Action action) {
  try {
    action();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

Record record_of(RecordType type, Metadata metadata = {}) {
  Record record;
  record.type = type;
  record.metadata = std::move(metadata);
  return record;
}

Record efficiency_status_of(const std::string &state) {
  return record_of(RecordType::efficiency_status,
                   Metadata{{"factors", std::string("atmosphere,qe")},
                            {"atmosphere.state", state},
                            {"atmosphere.source", std::string("test")},
                            {"qe.state", std::string("envelope")},
                            {"qe.source", std::string("test")}});
}

Record atmosphere_status_of(const std::string &state) {
  return record_of(RecordType::efficiency_status,
                   Metadata{{"factors", std::string("atmosphere")},
                            {"atmosphere.state", state},
                            {"atmosphere.source", std::string("test")}});
}

struct Photon {
  float x, y, cx, cy, t, w, zem, lambda;
};

Record photons_of(std::uint32_t telescope, const std::vector<Photon> &photons,
                  Frame frame = Frame::ground) {
  Record record;
  record.type = RecordType::photons;
  record.chunk.telescope_id = telescope;
  record.chunk.frame = frame;
  record.chunk.resize(photons.size(),
                      required_columns | column_bit(Column::zem) | column_bit(Column::wavelength));
  for (std::size_t i = 0; i < photons.size(); ++i) {
    record.chunk[Column::x][i] = photons[i].x;
    record.chunk[Column::y][i] = photons[i].y;
    record.chunk[Column::cx][i] = photons[i].cx;
    record.chunk[Column::cy][i] = photons[i].cy;
    record.chunk[Column::t][i] = photons[i].t;
    record.chunk[Column::w][i] = photons[i].w;
    record.chunk[Column::zem][i] = photons[i].zem;
    record.chunk[Column::wavelength][i] = photons[i].lambda;
  }
  return record;
}

Record telescope_of(std::int64_t id) {
  return record_of(RecordType::telescope,
                   Metadata{{"telescope_id", id}, {"array_index", std::int64_t{0}}});
}

Metadata run_metadata() {
  return Metadata{{"run_number", std::int64_t{5}},
                  {"observation_level_m", 2200.0},
                  {"energy_min_gev", 100.0},
                  {"energy_max_gev", 10000.0},
                  {"spectral_index", -2.0},
                  {"telescope_x_m", std::vector<double>{10.0, -30.0, 0.0}},
                  {"telescope_y_m", std::vector<double>{20.0, 40.0, 0.0}}};
}

Metadata event_metadata() {
  const double pi = std::numbers::pi;
  return Metadata{{"event_number", std::int64_t{1}},
                  {"primary_id", std::int64_t{14}},
                  {"energy_gev", 1000.0},
                  {"zenith_rad", pi / 2.0},
                  {"azimuth_rad", 0.0},
                  {"first_interaction_m", 25000.0},
                  {"core_x_m", std::vector<double>{100.0}},
                  {"core_y_m", std::vector<double>{50.0}}};
}

// Run with three telescopes: 1 and 2 carry one photon each (2 with weight `weight2`), 3 one.
std::vector<Record> stream_records(float weight2 = 1.0F, const std::string &atmosphere = "applied",
                                   float lambda = 400.0F) {
  std::vector<Record> records;
  records.push_back(record_of(RecordType::stream, Metadata{{"producer", std::string("test")}}));
  records.push_back(record_of(RecordType::run, run_metadata()));
  if (!atmosphere.empty()) {
    records.push_back(efficiency_status_of(atmosphere));
  }
  records.push_back(record_of(RecordType::event, event_metadata()));
  records.push_back(telescope_of(1));
  records.push_back(photons_of(1, {{1.5F, 0.5F, 0.0F, 0.5F, 12.5F, 1.0F, 20000.0F, lambda}}));
  records.push_back(
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{1}}}));
  records.push_back(telescope_of(2));
  records.push_back(photons_of(2, {{0.0F, 0.0F, 0.5F, 0.0F, 1.0F, weight2, 20000.0F, 500.0F}}));
  records.push_back(
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{2}}}));
  records.push_back(telescope_of(3));
  records.push_back(photons_of(3, {{0.0F, 0.0F, 0.0F, 0.0F, 1.0F, 1.0F, 20000.0F, 400.0F}}));
  records.push_back(
      record_of(RecordType::telescope_end, Metadata{{"telescope_id", std::int64_t{3}}}));
  records.push_back(record_of(RecordType::event_end));
  records.push_back(record_of(RecordType::run_end));
  return records;
}

std::string convert(const std::vector<Record> &records, GrisuOptions options,
                    std::string *warnings = nullptr) {
  std::ostringstream out;
  std::ostringstream warn;
  options.warnings = &warn;
  PsToGrisu converter{out, std::move(options)};
  for (const Record &record : records) {
    converter.write(record);
  }
  converter.finish();
  if (warnings != nullptr) {
    *warnings = warn.str();
  }
  return out.str();
}

std::vector<GrisuTelescope> tlloc_two() {
  std::istringstream in("# comment\n"
                        "* some free text\n"
                        "* TLLOC 7 -20.2 -10.1 0 0 0 0\n"
                        "* TLLOC 8 -40.0 30.0 5 0 0 0\n");
  return read_tlloc(in);
}

const char *const header_text = "* HEADF  <-- Start of header flag\n"
                                "* Produced by test-tool from a photon stream written by test\n"
                                "* CORSIKA run number: 5\n"
                                "* Energy range [TeV]: 0.1 10\n"
                                "* Spectral slope: -2\n"
                                "* Primary id: 14\n"
                                "* Primary zenith [deg]: 90\n"
                                "* Primary azimuth [deg, CORSIKA]: 0\n"
                                "* Observation height [m]: 2200\n"
                                "* DATAF  <-- end of header flag\n";

void transform_round_trips() {
  std::mt19937_64 rng{7};
  std::uniform_real_distribution<double> position(-1e4, 1e4);
  std::uniform_real_distribution<double> cosine(-0.6, 0.6);
  for (int i = 0; i < 1000; ++i) {
    const double x = position(rng);
    const double y = position(rng);
    const auto [gx, gy] = corsika_to_grisu_position(x, y);
    CHECK_NEAR(gx, -y, 0.0);
    CHECK_NEAR(gy, -x, 0.0);
    const auto [bx, by] = grisu_to_corsika_position(gx, gy);
    CHECK_NEAR(bx, x, 0.0);
    CHECK_NEAR(by, y, 0.0);
    const double cx = cosine(rng);
    const double cy = cosine(rng);
    const auto [dcx, dcy] = corsika_to_grisu_direction(cx, cy);
    CHECK_NEAR(dcx, -cy, 1e-12);
    CHECK_NEAR(dcy, -cx, 1e-12);
    const auto [rx, ry] = grisu_to_corsika_direction(dcx, dcy);
    CHECK_NEAR(rx, cx, 1e-12);
    CHECK_NEAR(ry, cy, 1e-12);
  }
  // Primary: CORSIKA azimuth 0 (towards +x) at zenith 90 deg is -y in GrISU.
  const auto [pcos, psin] = grisu_direction_from_angles(std::numbers::pi / 2.0, 0.0);
  CHECK_NEAR(pcos, 0.0, 1e-12);
  CHECK_NEAR(psin, -1.0, 1e-12);
  // The primary and photon rules agree: zenith 30 deg, azimuth 40 deg.
  const double ze = 30.0 * std::numbers::pi / 180.0;
  const double az = 40.0 * std::numbers::pi / 180.0;
  const auto [a, b] = grisu_direction_from_angles(ze, az);
  const auto [c, d] =
      corsika_to_grisu_direction(std::sin(ze) * std::cos(az), std::sin(ze) * std::sin(az));
  CHECK_NEAR(a, c, 1e-12);
  CHECK_NEAR(b, d, 1e-12);
}

void exact_text() {
  GrisuOptions options;
  options.producer = "test-tool";
  options.telescopes = tlloc_two();
  options.qe = 0.25;
  std::string warnings;
  const std::string text = convert(stream_records(), options, &warnings);
  const std::string expected = std::string(header_text) +
                               "R 0.25\n"
                               "H 2200\n"
                               "S 1 -50 -100 0 -1 25000 -1 -1 -1\n"
                               "P -20.5 -11.5 -0.5 +0 +20000 +12.5 400 3 7\n"
                               "P -40 +30 +0 -0.5 +20000 +1 500 3 8\n";
  CHECK(text == expected);
  if (text != expected) {
    std::fprintf(stderr, "got:\n%s\n", text.c_str());
  }
  CHECK(warnings.find("skipped 1 telescope") != std::string::npos);
}

void without_config_the_stream_id_is_used() {
  const std::string text = convert(stream_records(), GrisuOptions{});
  CHECK(text.find(" 3 1\nP ") != std::string::npos);
  CHECK(text.find(" 3 2\n") != std::string::npos);
  CHECK(text.find(" 3 3\n") != std::string::npos);
  CHECK(text.find("R 1\n") != std::string::npos);
}

void event_details_line() {
  GrisuOptions options;
  options.event_details = true;
  const std::string text = convert(stream_records(), options);
  CHECK(text.find("S 1 -50 -100 0 -1 25000 -1 -1 -1\nC 25000 -1 1\nP ") != std::string::npos);
}

void bunches() {
  GrisuOptions options;
  CHECK(throws_error([&] { convert(stream_records(3.0F), options); }));
  CHECK(throws_error([&] { convert(stream_records(0.5F), options); }));
  options.expand_bunches = true;
  const std::string text = convert(stream_records(3.0F), options);
  std::size_t copies = 0;
  for (std::size_t at = text.find(" 500 3 2\n"); at != std::string::npos;
       at = text.find(" 500 3 2\n", at + 1)) {
    ++copies;
  }
  CHECK(copies == 3);
  CHECK(throws_error([&] { convert(stream_records(0.5F), options); })); // still below 1
}

void tlloc_matching() {
  std::istringstream in("* TLLOC 1 -20.4 -10.0 0 0 0 0\n* TLLOC 2 -20.1 -10.0 0 0 0 0\n"
                        "* TLLOC 3 -40.6 30 0 0 0 0\n");
  GrisuOptions options;
  options.telescopes = read_tlloc(in);
  // Telescope 1 is at (-20, -10): ids 1 (0.4 m) and 2 (0.1 m) are inside 0.5 m, 2 is closer.
  // Telescope 2 is at (-40, 30): id 3 is 0.6 m away, so it is skipped.
  const std::string text = convert(stream_records(), options);
  CHECK(text.find("-11.5 -0.5 +0 +20000 +12.5 400 3 2\n") != std::string::npos);
  CHECK(text.find(" 500 3 ") == std::string::npos);
  std::istringstream twice("* TLLOC 1 0 0 0 0 0 0\n* TLLOC 1 1 1 0 0 0 0\n");
  CHECK(throws_error([&] { read_tlloc(twice); }));
  std::istringstream short_line("* TLLOC 1 0 0\n");
  CHECK(throws_error([&] { read_tlloc(short_line); }));
}

void error_cases() {
  GrisuOptions options;
  CHECK(throws_error([&] { convert(stream_records(1.0F, "pending"), options); }));
  CHECK(throws_error([&] { convert(stream_records(1.0F, "envelope"), options); }));
  CHECK(throws_error([&] { convert(stream_records(1.0F, ""), options); }));
  CHECK(throws_error([&] { convert(stream_records(1.0F, "applied", 0.0F), options); }));
  // Photons in the camera frame.
  std::vector<Record> camera = stream_records();
  camera[5] = photons_of(1, {{0.F, 0.F, 0.F, 0.F, 0.F, 1.F, 1.F, 400.F}}, Frame::camera);
  CHECK(throws_error([&] { convert(camera, options); }));
  // No observation height.
  std::vector<Record> no_height = stream_records();
  no_height[1].metadata.erase("observation_level_m");
  CHECK(throws_error([&] { convert(no_height, options); }));
  // No wavelength column.
  std::vector<Record> no_lambda = stream_records();
  no_lambda[5].chunk.resize(1, required_columns | column_bit(Column::zem));
  no_lambda[5].chunk[Column::w][0] = 1.0F;
  CHECK(throws_error([&] { convert(no_lambda, options); }));
  CHECK(throws_error([] {
    GrisuOptions bad;
    bad.qe = std::nan("");
    convert({}, bad);
  }));
  const float infinity = std::numeric_limits<float>::infinity();
  CHECK(throws_error([&] { convert(stream_records(1.0F, "applied", infinity), options); }));
  options.expand_bunches = true;
  CHECK(throws_error([&] { convert(stream_records(infinity), options); }));
  CHECK(throws_error([&] { convert(stream_records(1e9F), options); }));
  options.expand_bunches = false;
  // Telescope ids must be integers >= 1.
  for (const std::int64_t bad : {std::int64_t{0}, std::int64_t{-1}, std::int64_t{1} << 33}) {
    std::vector<Record> bad_id = stream_records();
    bad_id[4].metadata["telescope_id"] = bad;
    CHECK(throws_error([&] { convert(bad_id, options); }));
  }
  std::vector<Record> real_id = stream_records();
  real_id[4].metadata["telescope_id"] = 1.5;
  CHECK(throws_error([&] { convert(real_id, options); }));
  // A telescope efficiency status after the photons is rejected.
  std::vector<Record> late = stream_records();
  late.insert(late.begin() + 6, efficiency_status_of("pending"));
  CHECK(throws_error([&] { convert(late, options); }));
  // A telescope efficiency status overrides the run efficiency status.
  std::vector<Record> override_records = stream_records();
  override_records.insert(override_records.begin() + 5, efficiency_status_of("pending"));
  CHECK(throws_error([&] { convert(override_records, options); }));

  // Validate run and telescope transitions independently before applying the override.
  std::vector<Record> forward = stream_records(1.0F, "");
  forward.insert(forward.begin() + 2, atmosphere_status_of("applied"));
  forward.insert(forward.begin() + 2, atmosphere_status_of("envelope"));
  forward.insert(forward.begin() + 2, atmosphere_status_of("pending"));
  CHECK(!throws_error([&] { convert(forward, options); }));
  std::vector<Record> backward = stream_records(1.0F, "");
  backward.insert(backward.begin() + 2, atmosphere_status_of("pending"));
  backward.insert(backward.begin() + 2, atmosphere_status_of("envelope"));
  backward.insert(backward.begin() + 2, atmosphere_status_of("applied"));
  CHECK(throws_error([&] { convert(backward, options); }));
}

} // namespace

int main() {
  return run_checks([] {
    transform_round_trips();
    exact_text();
    without_config_the_stream_id_is_used();
    event_details_line();
    bunches();
    tlloc_matching();
    error_cases();
  });
}
