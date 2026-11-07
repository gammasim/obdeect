#include "check.hpp"

#include "dlamb/programs/octo.hpp"
#include "program_stream.hpp"

#include "gstrickd/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>

using namespace dlamb;

namespace {

// Geometry of octo.cc:107-124: 41 + 40 positions, times 0..20 and 20.5..40 ns.
void figure_eight_geometry() {
  const octo::Options o;
  const auto positions = octo::figure_eight(o);
  CHECK(positions.size() == 81);
  CHECK_NEAR(positions.front().time_ns, 0.0, 0.0);
  CHECK_NEAR(positions[40].time_ns, 20.0, 1e-12);
  CHECK_NEAR(positions[41].time_ns, 20.5, 1e-12);
  CHECK_NEAR(positions.back().time_ns, 40.0, 1e-12);
  for (std::size_t k = 0; k < positions.size(); ++k) {
    const Vec3 p = positions[k].position;
    CHECK_NEAR(p.z, 1000.0, 1e-9);
    // Each position lies on the circle of radius 20 m about (20, 0) or (-20, 0).
    const double cx = k <= 40 ? 20.0 : -20.0;
    CHECK_NEAR(std::hypot(p.x - cx, p.y), 20.0, 1e-9);
  }
  // Start of the first circle at the origin in plan view.
  CHECK_NEAR(positions.front().position.x, 0.0, 1e-9);
  CHECK_NEAR(positions.front().position.y, 0.0, 1e-9);
}

// Every stored photon reaches the sphere: its line passes within R of the centre.
// Direction is uniform in cos over the cone that contains the sphere (normalised).
void photons_hit_sphere_and_fill_cone() {
  const octo::Options o;
  const auto positions = octo::figure_eight(o);
  // The generation loop of octo::write, per source, to inspect the photons directly.
  std::uint64_t total = 0;
  double worst_miss = -1e300;
  double z_sum = 0.0, z_sq = 0.0;
  double expected_total = 0.0;
  for (std::size_t k = 0; k < positions.size(); ++k) {
    const Source s = octo::source_at(o, positions[k]);
    const std::uint64_t seed = source_seed(o.seed, k);
    const double cos_max = std::sqrt(1.0 - std::pow(12.5 / norm(positions[k].position), 2));
    const double u_max = 1.0 - cos_max;
    expected_total += expected_count(s);
    const std::uint64_t n = photon_count(s, seed, 1, 1);
    const Vec3 to_centre = normalised(-positions[k].position);
    for (std::uint64_t i = 0; i < n; ++i) {
      const Photon p = sample_photon(s, seed, 1, 1, i);
      ++total;
      // Distance from the centre to the ray: |P x d| for a unit direction d.
      worst_miss = std::max(worst_miss, norm(cross(p.position, p.direction)) - 12.5);
      const double z = (1.0 - dot(p.direction, to_centre)) / u_max;
      z_sum += z;
      z_sq += z * z;
      CHECK(p.weight == 1.0);
      CHECK_NEAR(p.time_ns, positions[k].time_ns, 1e-12);
      CHECK(p.wavelength_nm == 400.0);
    }
  }
  CHECK(worst_miss <= 1e-9);
  // Stored count against the expected count (Poisson-like).
  CHECK(std::abs(static_cast<double>(total) - expected_total) <= 4.0 * std::sqrt(expected_total));
  // z uniform on [0, 1]: mean 1/2 and variance 1/12.
  const double n = static_cast<double>(total);
  const double mean = z_sum / n;
  const double var = z_sq / n - mean * mean;
  CHECK(std::abs(mean - 0.5) <= 4.0 * std::sqrt(1.0 / 12.0 / n));
  CHECK(std::abs(var - 1.0 / 12.0) <= 4.0 * (1.0 / 12.0) * std::sqrt(2.0 / n));
}

// Whole stream: summary against the expected count, readable, sky frame, one event.
void stream_round_trip() {
  const octo::Options o;
  std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
  const StreamSummary summary = octo::write(stream, o);
  CHECK(summary.expected > 1000.0);
  CHECK(std::abs(static_cast<double>(summary.photons) - summary.expected) <=
        4.0 * std::sqrt(summary.expected));
  CHECK_NEAR(summary.weight, static_cast<double>(summary.photons), 1e-6);

  gstrickd::Reader reader(stream);
  gstrickd::Record record;
  int events = 0, telescopes = 0;
  std::uint64_t rows = 0;
  while (reader.next(record)) {
    if (record.type == gstrickd::RecordType::event)
      ++events;
    if (record.type == gstrickd::RecordType::telescope) {
      ++telescopes;
      CHECK(std::get<std::int64_t>(record.metadata.at("frame")) == 0);
      CHECK_NEAR(std::get<double>(record.metadata.at("r_m")), 12.5, 1e-12);
    }
    if (record.type == gstrickd::RecordType::photons) {
      CHECK(record.chunk.frame == gstrickd::Frame::sky);
      rows += record.chunk.rows();
    }
  }
  CHECK(reader.ended());
  CHECK(events == 1);
  CHECK(telescopes == 1);
  CHECK(rows == summary.photons);
}

// Emission altitude: every photon of the stream has zem = observation level + z (format.md).
void emission_altitude() {
  const octo::Options o;
  std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
  const StreamSummary summary = octo::write(stream, o);
  CHECK(summary.photons > 0);
  const StreamContents s = read_stream(stream.str());
  CHECK(!s.rows.empty());
  for (const auto &r : s.rows) {
    CHECK(r.has_zem);
    CHECK_NEAR(r.zem - r.z, o.observation_level_m, 1e-2);
  }
}

void validation() {
  octo::Options bad;
  bad.distance_m = 10.0; // inside the sphere
  bool rejected = false;
  try {
    octo::validate(bad);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
}

} // namespace

int main() {
  return run_checks([] {
    figure_eight_geometry();
    photons_hit_sphere_and_fill_cone();
    stream_round_trip();
    emission_altitude();
    validation();
  });
}
