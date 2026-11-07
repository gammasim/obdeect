#include "check.hpp"
#include "program_stream.hpp"

#include "dlamb/programs/fpls.hpp"

#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>

using namespace dlamb;
using namespace dlamb::programs;

namespace {

std::string stream_bytes(const fpls::Options &options) {
  std::ostringstream out;
  fpls::write_stream(out, options);
  return out.str();
}

void rotation_and_steps() {
  // Parallel beam (no angular spread): each time step t = k * 0.25 ns points along the
  // LED axis rotated by 2 pi t / period about z. Period 2 ns gives 9 steps (t = 0 ... 2).
  fpls::Options options;
  options.rotation_period_ns = 2.0;
  options.photons = 1000.0;
  options.angular = "Parallel";
  const auto s = read_stream(stream_bytes(options));

  CHECK(count_type(s, gstrickd::RecordType::event) == 1);
  CHECK(s.rows.size() == 9 * 1000);
  // Axis: DirectionVect(-75 deg, 0) = (cos 75 deg ... ): cx = cos(-75 deg), cz = sin(-75 deg).
  const double ax = std::cos(-75.0 * pi / 180.0);
  const double az = std::sin(-75.0 * pi / 180.0);
  for (const Row &r : s.rows) {
    CHECK(r.frame == 0); // sky
    CHECK_NEAR(r.x, 0.0, 1e-6);
    CHECK_NEAR(r.y, 0.0, 1e-6);
    CHECK_NEAR(r.z, 0.7 * 1.4, 1e-6); // LED height 0.7 R
    CHECK_NEAR(r.zem, 2100.0 + 0.98, 1e-2);
    const double k = std::round(static_cast<double>(r.t) / 0.25);
    CHECK_NEAR(r.t, 0.25 * k, 1e-5);
    const double phi = 2.0 * pi * (0.25 * k) / 2.0;
    CHECK_NEAR(r.cx, ax * std::cos(phi), 1e-5);
    CHECK_NEAR(r.cy, ax * std::sin(phi), 1e-5);
    CHECK_NEAR(r.cz, az, 1e-5);
    CHECK(r.w == 1.0F);
  }
}

void isotropic_counts() {
  // Default beam (isotropic): 5000 photons per step; period 0.5 ns gives steps 0, 0.25, 0.5.
  fpls::Options options;
  options.rotation_period_ns = 0.5;
  const auto s = read_stream(stream_bytes(options));
  CHECK(s.rows.size() == 3 * 5000);
  for (const Row &r : s.rows) {
    CHECK_NEAR(r.cx * r.cx + r.cy * r.cy + r.cz * r.cz, 1.0, 1e-5);
    CHECK(r.wavelength == 400.0F);
  }
}

} // namespace

int main() {
  return run_checks([] {
    rotation_and_steps();
    isotropic_counts();
  });
}
