#include "check.hpp"
#include "sky_stream_check.hpp"

#include "dlamb/programs/xyls.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <variant>

using namespace dlamb;
using namespace dlamb::programs;
using gstrickd::Column;

namespace {

// Number of scan positions of xyls.cc for a camera radius r (cm): the loops of the original.
std::uint64_t scan_positions(double x0, double y0, double r) {
  std::uint64_t n = 0;
  for (double x = x0; x < x0 + 0.4 * r; x += 0.001 * r)
    ++n;
  for (double y = y0; y < y0 + 0.2 * r; y += 0.001 * r)
    ++n;
  return n;
}

// Default: 100 photons per position, 1 bunch of size 1, parallel beam at height 6.15 cm.
void default_scan() {
  const XylsOptions o;
  const auto plan = plan_xyls(o);
  const std::uint64_t positions = scan_positions(o.x_cm, o.y_cm, o.radius_cm);
  CHECK(positions > 500 && positions < 700);
  CHECK(plan.cases.size() == 1 && plan.cases[0].size() == positions);
  // Defaults derived from 123 cm: x = 36.9 cm, y = 24.6 cm, height 6.15 cm.
  CHECK_NEAR(o.x_cm, 36.9, 1e-12);
  CHECK_NEAR(o.height_cm, 6.15, 1e-12);

  const auto s = read_stream(write_to_string([&](std::ostream &out) { write_xyls(out, o); }));
  CHECK(s.ended && s.events == 1);
  CHECK(s.photons.at(1) == 100 * positions);
  for (const auto &c : s.chunks) {
    CHECK(c.frame == gstrickd::Frame::sky);
    for (std::size_t i = 0; i < c.rows(); ++i) {
      CHECK_NEAR(c[Column::cx][i], 0.0, 1e-6);
      CHECK_NEAR(c[Column::cy][i], 0.0, 1e-6);
      CHECK_NEAR(c[Column::cz][i], -1.0, 1e-6);
      CHECK_NEAR(c[Column::w][i], 1.0, 1e-6);
      CHECK_NEAR(c[Column::z][i], 0.0615, 1e-6);
    }
  }
  // The first photon sits at the scan start.
  const auto &first = std::get<Emitter>(plan.cases[0][0].source.geometry).position;
  CHECK_NEAR(first.x, 0.369, 1e-12);
  CHECK_NEAR(first.y, 0.246, 1e-12);
}

// Non-default: radius 200 cm, height 10 cm, 3 photons per bunch of size 2 per position.
void non_default_scan() {
  XylsOptions o;
  o.radius_cm = 200.0;
  o.x_cm = 10.0;
  o.y_cm = 5.0;
  o.height_cm = 10.0;
  o.bunch_size = 2.0;
  o.photons_per_position = 6.0;
  const std::uint64_t positions = scan_positions(o.x_cm, o.y_cm, o.radius_cm);
  const auto s = read_stream(write_to_string([&](std::ostream &out) { write_xyls(out, o); }));
  CHECK(s.photons.at(1) == 3 * positions);
  for (const auto &c : s.chunks)
    for (std::size_t i = 0; i < c.rows(); ++i) {
      CHECK_NEAR(c[Column::w][i], 2.0, 1e-6);
      CHECK_NEAR(c[Column::z][i], 0.10, 1e-6);
      CHECK(c[Column::x][i] >= 0.10 - 1e-6 && c[Column::x][i] < 0.10 + 0.8 + 1e-6);
      CHECK(c[Column::y][i] >= 0.05 - 1e-6 && c[Column::y][i] < 0.05 + 0.4 + 1e-6);
    }
}

// Negative photon counts per position are rejected, not clamped to zero.
void rejects_negative_photons() {
  XylsOptions o;
  o.photons_per_position = -1.0;
  bool threw = false;
  try {
    (void)plan_xyls(o);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw);
}

// A zero, negative, NaN or infinite radius is rejected (it used to give an empty stream).
void rejects_bad_radius() {
  const double radii[] = {0.0, -1.0, -123.0, std::nan(""), std::numeric_limits<double>::infinity()};
  for (const double r : radii) {
    XylsOptions o;
    o.radius_cm = r;
    bool threw = false;
    try {
      (void)plan_xyls(o);
    } catch (const std::invalid_argument &) {
      threw = true;
    }
    CHECK(threw);
  }
  XylsOptions ok;
  ok.radius_cm = 1e-3; // tiny but positive: accepted, with a non-empty scan
  CHECK(!plan_xyls(ok).cases[0].empty());
}

} // namespace

int main() {
  return run_checks([] {
    default_scan();
    non_default_scan();
    rejects_negative_photons();
    rejects_bad_radius();
  });
}
