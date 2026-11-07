#include "check.hpp"
#include "program_stream.hpp"

#include "dlamb/programs/ff-1m.hpp"

#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

using namespace dlamb;
using namespace dlamb::programs;

namespace {

std::string stream_bytes(const ff_1m::Options &options) {
  std::ostringstream out;
  ff_1m::write_stream(out, options);
  return out.str();
}

// u = 1 - cos(theta) of the fiducial cone of LightSource::Emit: sin(theta) = R / d
// (measured on the original ff-1m output, scripts/compare_ff-1m.sh).
double cone_umax(double distance_m, double radius_m) {
  const double sine = radius_m / distance_m;
  return 1.0 - std::sqrt(1.0 - sine * sine);
}

void default_options() {
  ff_1m::Options options;
  options.events = {20};
  const auto s = read_stream(stream_bytes(options));

  // Structure: one run, 20 events each with telescope, photons, telescope_end, event_end.
  CHECK(s.types.size() > 2);
  CHECK(count_type(s, gstrickd::RecordType::run) == 1);
  CHECK(count_type(s, gstrickd::RecordType::event) == 20);
  CHECK(count_type(s, gstrickd::RecordType::telescope) == 20);
  CHECK(count_type(s, gstrickd::RecordType::telescope_end) == 20);
  CHECK(count_type(s, gstrickd::RecordType::event_end) == 20);
  CHECK(count_type(s, gstrickd::RecordType::run_end) == 1);

  // Count: 2.5e6 photons into 4 pi, of which the fiducial cone keeps
  // N = 2.5e6 * (1 - cos theta) / 2 per event.
  const double distance = 16.0;
  const double umax = cone_umax(distance, 1.2);
  const double expected = 2.5e6 * umax / 2.0;
  CHECK(s.rows.size() > 0);
  CHECK_NEAR(static_cast<double>(s.rows.size()), 20.0 * expected, 4.0 * std::sqrt(20.0 * expected));

  // Frame sky, unit weights, all at the same time and wavelength, emitted from the
  // flasher at (0, 0, 16 m); zem = observation level 2150 m + 16 m.
  // Direction: the cone points at the telescope centre, so u = 1 + cz is uniform on
  // [0, umax] (solid-angle uniform), mean umax / 2, variance umax^2 / 12.
  double sum = 0.0;
  double sum2 = 0.0;
  for (const Row &r : s.rows) {
    CHECK(r.frame == 0);
    CHECK(r.w == 1.0F);
    CHECK(r.t == 0.0F);
    CHECK(r.wavelength == 400.0F);
    CHECK(r.x == 0.0F && r.y == 0.0F);
    CHECK_NEAR(r.z, 16.0, 1e-4);
    CHECK(r.has_zem);
    CHECK_NEAR(r.zem, 2166.0, 1e-2);
    CHECK_NEAR(r.cx * r.cx + r.cy * r.cy + r.cz * r.cz, 1.0, 1e-5);
    const double u = 1.0 + r.cz;
    CHECK(u >= -1e-5 && u <= umax + 1e-5);
    sum += u;
    sum2 += u * u;
  }
  const double n = static_cast<double>(s.rows.size());
  const double mean = sum / n;
  const double variance = sum2 / n - mean * mean;
  CHECK_NEAR(mean, umax / 2.0, 4.0 * std::sqrt(umax * umax / 12.0 / n));
  CHECK_NEAR(variance, umax * umax / 12.0, 4.0 * umax * umax / 12.0 * std::sqrt(2.0 / n));
}

void non_default_options() {
  // Offset flasher (1 m, 0 m, 10 m), bunches of 2 photons, Simple:4 pulse (time spread
  // +-2 ns), Cherenkov 300-600 nm spectrum, 1e6 photons into 4 pi per event.
  ff_1m::Options options;
  options.x_cm = 100.0;
  options.distance_cm = 1000.0;
  options.bunch_size = 2.0;
  options.lightpulse = "Simple:4";
  options.spectrum = "Cherenkov:300:600";
  options.photons = {1e6};
  options.events = {1};
  options.altitude_m = 2000.0;
  const auto s = read_stream(stream_bytes(options));

  const double distance = std::sqrt(1.0 + 100.0);
  const double umax = cone_umax(distance, 1.2);
  const double emitted = 1e6 / 2.0 * umax / 2.0; // bunches of 2 photons
  // Weight sum = 2 photons per bunch: mean 2 * emitted, Poisson-like count.
  double weight = 0.0;
  double time_sum = 0.0;
  double time_sum2 = 0.0;
  double wl_sum = 0.0;
  double wl_sum2 = 0.0;
  double u_sum = 0.0;
  const Vec3 axis{-1.0 / distance, 0.0, -10.0 / distance}; // towards the telescope centre
  for (const Row &r : s.rows) {
    weight += r.w;
    CHECK_NEAR(r.x, 1.0, 1e-4);
    CHECK_NEAR(r.y, 0.0, 1e-4);
    CHECK_NEAR(r.z, 10.0, 1e-4);
    CHECK_NEAR(r.zem, 2010.0, 1e-2);
    CHECK(r.w == 2.0F);
    time_sum += r.t;
    time_sum2 += static_cast<double>(r.t) * r.t;
    wl_sum += r.wavelength;
    wl_sum2 += static_cast<double>(r.wavelength) * r.wavelength;
    const double u = 1.0 - (r.cx * axis.x + r.cy * axis.y + r.cz * axis.z);
    CHECK(u >= -1e-5 && u <= umax + 1e-5);
    u_sum += u;
  }
  const double n = static_cast<double>(s.rows.size());
  CHECK_NEAR(weight, 2.0 * emitted, 4.0 * 2.0 * std::sqrt(emitted));
  CHECK_NEAR(u_sum / n, umax / 2.0, 4.0 * std::sqrt(umax * umax / 12.0 / n));
  // Simple:4 is a top hat of full width 4 ns centred on 0 (variance 16/12).
  const double t_mean = time_sum / n;
  CHECK_NEAR(t_mean, 0.0, 4.0 * std::sqrt(16.0 / 12.0 / n));
  // Cherenkov spectrum: density proportional to 1/lambda^2 on [300, 600] nm.
  const InverseSquare spectrum{300.0, 600.0};
  const double wl_mean = wl_sum / n;
  const double wl_var = wl_sum2 / n - wl_mean * wl_mean;
  CHECK_NEAR(wl_mean, spectrum.mean(), 4.0 * std::sqrt(spectrum.variance() / n));
  CHECK_NEAR(wl_var, spectrum.variance(), 4.0 * spectrum.variance() * std::sqrt(2.0 / n));
  for (const Row &r : s.rows) {
    CHECK(r.wavelength >= 300.0F - 1e-3F && r.wavelength <= 600.0F + 1e-3F);
  }
  (void)time_sum2;
}

void event_numbering_and_repeat() {
  // Two photon levels: events {1, 3}; the second level has 3 events (events[1]).
  ff_1m::Options options;
  options.photons = {1e5, 1e5};
  options.events = {1, 3};
  const auto s = read_stream(stream_bytes(options));
  CHECK(count_type(s, gstrickd::RecordType::event) == 4);
  std::int64_t last_event = 0;
  for (const auto &[type, metadata] : s.metadata) {
    if (type == gstrickd::RecordType::event) {
      const auto number = std::get<std::int64_t>(metadata.at("event_number"));
      CHECK(number == last_event + 1);
      last_event = number;
    }
  }
  CHECK(last_event == 4);
}

} // namespace

int main() {
  return run_checks([] {
    default_options();
    non_default_options();
    event_numbering_and_repeat();
  });
}
