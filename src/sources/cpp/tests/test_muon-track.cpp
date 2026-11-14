#include "check.hpp"

#include "obdeect/sources/emission_stream.hpp"
#include "obdeect/sources/programs/muon-track.hpp"
#include "obdeect/sources/track.hpp"
#include "program_stream.hpp"

#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace obdeect::sources;

namespace {

// Stored photons of all tracks of one event.
std::vector<Photon> stored_photons(const std::vector<Track> &tracks, std::uint64_t seed,
                                   std::uint32_t event = 1) {
  std::vector<Photon> out;
  for (std::size_t k = 0; k < tracks.size(); ++k)
    emit_track(tracks[k], source_seed(seed, k), event, 1, out);
  return out;
}

double sum_weight(const std::vector<Photon> &p) {
  double s = 0.0;
  for (const auto &x : p)
    s += x.weight;
  return s;
}

double sum_weight_sq(const std::vector<Photon> &p) {
  double s = 0.0;
  for (const auto &x : p)
    s += x.weight * x.weight;
  return s;
}

// True when f throws std::invalid_argument.
template <class F> bool throws_invalid(F f) {
  try {
    f();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

// Track of the muon-track default geometry but with a given sphere radius and index.
Track vertical_track(double length, double cos_ct, double radius) {
  Track t;
  t.start = {0.0, 0.0, length};
  t.direction = {0.0, 0.0, -1.0};
  t.length_m = length;
  t.cos_theta = cos_ct;
  t.spectrum = InverseSquare{300.0, 600.0};
  t.photons = 1e5;
  t.target = SphereTarget{{0.0, 0.0, 0.0}, radius};
  return t;
}

// Arc kernel against brute force over the azimuth: the fraction of azimuths that hit.
void arc_matches_brute_force() {
  Track t = vertical_track(1000.0, std::cos(1.0 * pi / 180.0), 15.0);
  t.start = {20.0, -5.0, 1000.0};
  t.direction = normalised(Vec3{-0.2, 0.1, -1.0});
  t.cos_theta = std::cos(19.0 * pi / 180.0); // the sphere lies at 19 deg: partial arc
  const double s = 400.0;
  const TrackArc arc = track_arc(t, s);
  const Vec3 point = t.start + t.direction * s;
  const Basis b = orthonormal_basis(t.direction);
  const double sin_t = std::sin(std::acos(t.cos_theta));
  const std::size_t n = 200000;
  std::size_t hits = 0;
  for (std::size_t k = 0; k < n; ++k) {
    const double phi = two_pi * (static_cast<double>(k) + 0.5) / n;
    const Vec3 e =
        t.direction * t.cos_theta + (b.e1 * std::cos(phi) + b.e2 * std::sin(phi)) * sin_t;
    const Vec3 to = normalised(t.target.centre - point);
    if (dot(e, to) > std::sqrt(1.0 - std::pow(15.0 / norm(t.target.centre - point), 2)))
      ++hits;
  }
  const double fraction = static_cast<double>(hits) / n;
  CHECK(arc.hit);
  CHECK_NEAR(arc.half / pi, fraction, 2e-5);
}

// Axis point: all azimuths hit if the emission angle is inside the sphere's cone.
void arc_on_axis() {
  Track t = vertical_track(3000.0, std::cos(0.2 * pi / 180.0), 15.0); // inside asin(15/3000)
  CHECK(track_arc(t, 0.0).hit);
  CHECK_NEAR(track_arc(t, 0.0).half, pi, 1e-12);
  t.cos_theta = std::cos(1.0 * pi / 180.0); // wider than asin(15/3000) = 0.286 deg
  CHECK(!track_arc(t, 0.0).hit);
  // A point inside the sphere is hit in every direction.
  Track inside = vertical_track(3000.0, std::cos(60.0 * pi / 180.0), 15.0);
  CHECK(track_arc(inside, 2995.0).hit);
  CHECK_NEAR(track_arc(inside, 2995.0).half, pi, 1e-12);
}

// Default options: Frank-Tamm yield, counts and weight sums against the analytic values.
void default_options() {
  muon_track::Options o;
  const auto segments = muon_track::tracks(o);
  CHECK(segments.size() == 1);
  const double sin2 = 1.0 - std::pow(1.0 / (1.0 + o.index_minus_one), 2);
  const double per_m =
      2.0 * pi * muon_track::alpha_fine * sin2 * (1.0 / 240.0 - 1.0 / 1000.0) * 1e9;
  CHECK_NEAR(segments[0].photons, per_m * 3000.0, 1e-9 * per_m * 3000.0);
  CHECK_NEAR(segments[0].time_ns, -3000.0 / speed_of_light_m_per_ns, 1e-9);
  // Emitted count (4 pi) is Poisson-like around the mean.
  const std::uint64_t emitted = track_count(segments[0], 1, 1, 1);
  CHECK(std::abs(static_cast<double>(emitted) - segments[0].photons) <= 1.0);
  // Expected weight of a central track (impact at the telescope) is positive and matches
  // the stored sum within 4 sigma.
  const double expected = expected_track_weight(segments[0]);
  CHECK(expected > 0.0);
  const auto photons = stored_photons(segments, 7);
  const double sigma = std::sqrt(sum_weight_sq(photons));
  CHECK(std::abs(sum_weight(photons) - expected) <= 4.0 * sigma);
  for (const auto &p : photons) {
    CHECK(p.weight >= 0.0 && p.weight <= 1.0);
    CHECK(p.time_ns <= 0.0);
  }
}

// Non-default options: angles, impact, index, band, altitude and radius. The emission angle
// is exact, the wavelength follows 1/lambda^2, and the weight sum matches the quadrature.
void non_default_options() {
  muon_track::Options o;
  o.theta_deg = 20.0;
  o.phi_deg = 30.0;
  o.impact_x_m = 40.0;
  o.impact_y_m = -30.0;
  o.h_start_m = 2000.0;
  o.wavelength_min_nm = 300.0;
  o.wavelength_max_nm = 600.0;
  o.telescope_radius_m = 15.0;
  o.index_minus_one = 0.00028;
  const auto segments = muon_track::tracks(o);
  CHECK(segments.size() == 1);
  const double dist = 2000.0 / std::cos(20.0 * pi / 180.0);
  CHECK_NEAR(segments[0].length_m, dist, 1e-9);
  // The track ends at the impact point on the observation level.
  CHECK_NEAR(segments[0].start.x + segments[0].direction.x * dist, 40.0, 1e-7);
  CHECK_NEAR(segments[0].start.y + segments[0].direction.y * dist, -30.0, 1e-7);

  const double cos_ct = 1.0 / (1.0 + o.index_minus_one);
  const double expected = expected_track_weight(segments[0]);
  const auto photons = stored_photons(segments, 11);
  CHECK(photons.size() > 200);
  const double sigma = std::sqrt(sum_weight_sq(photons));
  CHECK(std::abs(sum_weight(photons) - expected) <= 4.0 * sigma);

  // Emission angle to the track axis is exactly theta_c.
  double worst = 0.0;
  for (const auto &p : photons)
    worst = std::max(worst, std::abs(dot(p.direction, segments[0].direction) - cos_ct));
  CHECK(worst < 1e-12);

  // Wavelength: 1/lambda^2 moments within 4 sigma.
  const InverseSquare law{300.0, 600.0};
  double s1 = 0.0, s2 = 0.0;
  for (const auto &p : photons) {
    s1 += p.wavelength_nm;
    s2 += p.wavelength_nm * p.wavelength_nm;
  }
  const double n = static_cast<double>(photons.size());
  const double mean = s1 / n;
  const double var = s2 / n - mean * mean;
  CHECK(std::abs(mean - law.mean()) <= 4.0 * std::sqrt(law.variance() / n));
  CHECK(std::abs(var - law.variance()) <= 4.0 * law.variance() * std::sqrt(2.0 / n));
  for (const auto &p : photons)
    CHECK(p.wavelength_nm >= 300.0 && p.wavelength_nm <= 600.0);

  // Emission times lie between the start (-dist/c) and the arrival at the observation level.
  for (const auto &p : photons)
    CHECK(p.time_ns >= -dist / speed_of_light_m_per_ns - 1e-9 && p.time_ns <= 1e-9);
}

// Profile option: the track is cut into 10 m steps with the index at each mid-point, and the
// steps cover the whole track.
void profile_steps() {
  muon_track::Options o;
  o.h_start_m = 3000.0;
  o.index_profile.altitude_m = {0.0, 3000.0};
  o.index_profile.index_minus_one = {0.000292, 0.000220};
  const auto segments = muon_track::tracks(o);
  CHECK(segments.size() == 300);
  double length = 0.0;
  for (const auto &t : segments)
    length += t.length_m;
  CHECK_NEAR(length, 3000.0, 1e-7);
  // The first step (lowest, at the observation level) uses the index at 5 m altitude.
  const double n0 = 1.0 + 0.000292 + (0.000220 - 0.000292) * 5.0 / 3000.0;
  CHECK_NEAR(segments.front().cos_theta, 1.0 / n0, 1e-12);
  CHECK_NEAR(segments.front().time_ns, -10.0 / speed_of_light_m_per_ns, 1e-9);
  // Rejected: a profile with unsorted altitudes.
  muon_track::Options bad = o;
  bad.index_profile.altitude_m = {3000.0, 0.0};
  bool rejected = false;
  try {
    muon_track::validate(bad);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
}

// Options that are not supported are rejected.
void rejects_unsupported() {
  muon_track::Options o;
  o.multiple_scattering = true;
  bool rejected = false;
  try {
    muon_track::validate(o);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
  o = muon_track::Options{};
  o.rscat_m = 1.0;
  rejected = false;
  try {
    muon_track::validate(o);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
}

// Stream: readable by obdeect::stream::Reader; run, events, sky frame and totals agree.
void stream_round_trip() {
  muon_track::Options o;
  o.theta_deg = 20.0;
  o.impact_x_m = 40.0;
  o.impact_y_m = -30.0;
  o.h_start_m = 2000.0;
  o.index_minus_one = 0.00028;
  o.events = 2;
  std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
  const StreamSummary summary = muon_track::write(stream, o);
  CHECK(summary.photons > 0);
  CHECK_NEAR(summary.weight, summary.expected, 4.0 * std::sqrt(summary.expected + 1.0));

  obdeect::stream::Reader reader(stream);
  obdeect::stream::Record record;
  int runs = 0, events = 0, telescopes = 0;
  std::uint64_t rows = 0;
  double telescope_weight = 0.0;
  while (reader.next(record)) {
    switch (record.type) {
    case obdeect::stream::RecordType::run:
      ++runs;
      CHECK(std::get<std::int64_t>(record.metadata.at("run_number")) == 1);
      CHECK(std::get<std::int64_t>(record.metadata.at("refraction_applied")) == 0);
      CHECK(std::get<std::int64_t>(record.metadata.at("curvature_applied")) == 0);
      break;
    case obdeect::stream::RecordType::event:
      ++events;
      break;
    case obdeect::stream::RecordType::telescope:
      ++telescopes;
      CHECK(std::get<std::int64_t>(record.metadata.at("frame")) ==
            static_cast<std::int64_t>(obdeect::stream::Frame::sky));
      CHECK_NEAR(std::get<double>(record.metadata.at("r_m")), 15.0, 1e-12);
      break;
    case obdeect::stream::RecordType::photons:
      CHECK(record.chunk.frame == obdeect::stream::Frame::sky);
      CHECK(record.chunk.rows() <= obdeect::stream::max_chunk_rows);
      rows += record.chunk.rows();
      break;
    case obdeect::stream::RecordType::telescope_end:
      telescope_weight += std::get<double>(record.metadata.at("photons"));
      break;
    default:
      break;
    }
  }
  CHECK(reader.ended());
  CHECK(runs == 1);
  CHECK(events == 2);
  CHECK(telescopes == 2);
  CHECK(rows == summary.photons);
  CHECK_NEAR(telescope_weight, summary.weight, 1e-3 * summary.weight + 1e-6);
}

// Emission altitude: every photon of the stream has zem = observation level + z (format.md).
void emission_altitude() {
  muon_track::Options o;
  o.theta_deg = 20.0;
  o.impact_x_m = 40.0;
  o.impact_y_m = -30.0;
  o.h_start_m = 2000.0;
  o.index_minus_one = 0.00028;
  o.events = 1;
  std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
  const StreamSummary summary = muon_track::write(stream, o);
  CHECK(summary.photons > 0);
  const StreamContents s = read_stream(stream.str());
  CHECK(!s.rows.empty());
  for (const auto &r : s.rows) {
    CHECK(r.has_zem);
    CHECK_NEAR(r.zem - r.z, o.observation_level_m, 1e-2);
  }
}

// Track magnitudes beyond 1e30 are rejected, as in Source::validate (finite float32 columns).
void rejects_large_magnitudes() {
  const double cos_ct = std::cos(1.0 * pi / 180.0);
  CHECK(!throws_invalid([&] { validate(vertical_track(1000.0, cos_ct, 15.0)); }));
  Track t = vertical_track(1000.0, cos_ct, 15.0);
  t.start.x = 1e100;
  CHECK(throws_invalid([&] { validate(t); }));
  t = vertical_track(1000.0, cos_ct, 15.0);
  t.length_m = 1e100;
  CHECK(throws_invalid([&] { validate(t); }));
  t = vertical_track(1000.0, cos_ct, 15.0);
  t.target.centre.y = 1e100;
  CHECK(throws_invalid([&] { validate(t); }));
  t = vertical_track(1000.0, cos_ct, 1e100);
  CHECK(throws_invalid([&] { validate(t); }));
  t = vertical_track(1000.0, cos_ct, 15.0);
  t.time_ns = -1e100;
  CHECK(throws_invalid([&] { validate(t); }));
  t = vertical_track(1000.0, cos_ct, 15.0);
  t.spectrum = Fixed{1e100};
  CHECK(throws_invalid([&] { validate(t); }));
  t.spectrum = InverseSquare{300.0, 1e100};
  CHECK(throws_invalid([&] { validate(t); }));
}

// Track photon counts above 2^53 are rejected (track_count would convert them to uint64_t
// without loss); 2^53 itself is accepted, as in Source::validate.
void rejects_large_photon_counts() {
  Track t = vertical_track(1000.0, std::cos(1.0 * pi / 180.0), 15.0);
  t.photons = max_photon_count;
  CHECK(!throws_invalid([&] { validate(t); }));
  t.photons = 1e16;
  CHECK(throws_invalid([&] { validate(t); }));
  std::vector<Photon> out;
  CHECK(throws_invalid([&] { emit_track(t, 1, 1, 1, out); }));
  t.photons = std::nan("");
  CHECK(throws_invalid([&] { validate(t); }));
}

// CLI value checks of muon-track: event count range, primary names and index profile file.
void cli_values() {
  CHECK(throws_invalid([] { (void)muon_track::event_count(-1); }));
  CHECK(throws_invalid([] { (void)muon_track::event_count(0); }));
  CHECK(throws_invalid([] { (void)muon_track::event_count(4294967296LL); }));
  CHECK(muon_track::event_count(1) == 1);
  CHECK(muon_track::event_count(4294967295LL) == 4294967295U);

  // mu+ is CORSIKA 5 (charge +1), mu- is 6 (charge -1); "muon" is mu-.
  double mass = 0.0, charge = 0.0;
  CHECK(muon_track::primary_id_from_name("mu+") == 5);
  muon_track::mass_charge(muon_track::primary_id_from_name("mu+"), mass, charge);
  CHECK(charge == 1.0);
  CHECK(muon_track::primary_id_from_name("mu-") == 6);
  CHECK(muon_track::primary_id_from_name("muon") == 6);
  muon_track::mass_charge(muon_track::primary_id_from_name("muon"), mass, charge);
  CHECK(charge == -1.0);
  CHECK(throws_invalid([] { (void)muon_track::primary_id_from_name("tau"); }));

  // An empty index profile is an error, not the constant index.
  std::istringstream empty("");
  CHECK(throws_invalid([&] { (void)muon_track::read_index_profile(empty, "empty.txt"); }));
  std::istringstream garbage("0 0.000292\nabc\n");
  CHECK(throws_invalid([&] { (void)muon_track::read_index_profile(garbage, "garbage.txt"); }));
  std::istringstream dangling("0 0.000292\n3000\n");
  CHECK(throws_invalid([&] { (void)muon_track::read_index_profile(dangling, "dangling.txt"); }));
  std::istringstream two("0 0.000292\n3000 0.000220\n");
  const muon_track::IndexProfile profile = muon_track::read_index_profile(two, "two.txt");
  CHECK(profile.altitude_m.size() == 2);
  CHECK_NEAR(profile.index_minus_one[1], 0.000220, 1e-15);
}

// The default primary id, charge and mass describe the same particle (mu-, CORSIKA 6).
void default_primary_matches_charge() {
  const muon_track::Options o;
  double mass = 0.0, charge = 0.0;
  muon_track::mass_charge(o.primary_id, mass, charge);
  CHECK(o.primary_id == 6);
  CHECK(charge == o.charge);
  CHECK(mass == o.mass_gev);
  CHECK(muon_track::primary_id_from_name("muon") == o.primary_id);
}

// A profile with index values but no altitudes is rejected, not run as constant index.
void rejects_profile_size_mismatch() {
  muon_track::Options o;
  o.index_profile.index_minus_one = {0.000292};
  CHECK(throws_invalid([&] { muon_track::validate(o); }));
  CHECK(throws_invalid([&] { (void)muon_track::tracks(o); }));
  o = muon_track::Options{};
  o.index_profile.altitude_m = {0.0, 3000.0};
  CHECK(throws_invalid([&] { muon_track::validate(o); }));
  o = muon_track::Options{};
  CHECK(!throws_invalid([&] { muon_track::validate(o); }));
}

// --momentum gives the same beta whether it comes before or after --primary and -m.
void momentum_option_order() {
  const auto before =
      muon_track::parse_command_line({"--momentum", "1", "--primary", "electron"}).options;
  const auto after =
      muon_track::parse_command_line({"--primary", "electron", "--momentum", "1"}).options;
  const auto mass_last =
      muon_track::parse_command_line({"-p", "1", "--primary", "electron", "-m", "0.000511"})
          .options;
  const double expected = 1.0 / std::sqrt(1.0 + 0.000511 * 0.000511);
  CHECK(before.mass_gev == 0.000511);
  CHECK_NEAR(before.beta, expected, 1e-15);
  CHECK_NEAR(after.beta, before.beta, 0.0);
  CHECK_NEAR(mass_last.beta, before.beta, 0.0);
  // --beta and --momentum both set the speed: rejected in either order.
  CHECK(throws_invalid(
      [] { (void)muon_track::parse_command_line({"--beta", "0.5", "--momentum", "1"}); }));
  CHECK(throws_invalid(
      [] { (void)muon_track::parse_command_line({"--momentum", "1", "--beta", "0.5"}); }));
  CHECK(throws_invalid([] { (void)muon_track::parse_command_line({"-p", "0"}); }));
}

} // namespace

int main() {
  return run_checks([] {
    arc_matches_brute_force();
    arc_on_axis();
    default_options();
    non_default_options();
    profile_steps();
    rejects_unsupported();
    stream_round_trip();
    emission_altitude();
    rejects_large_magnitudes();
    rejects_large_photon_counts();
    cli_values();
    default_primary_matches_charge();
    rejects_profile_size_mismatch();
    momentum_option_order();
  });
}
