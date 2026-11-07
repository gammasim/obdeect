#include "check.hpp"
#include "program_stream.hpp"

#include "dlamb/programs/ff-gct.hpp"

#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>

using namespace dlamb;
using namespace dlamb::programs;

namespace {

std::string stream_bytes(const ff_gct::Options &options) {
  std::ostringstream out;
  ff_gct::write_stream(out, options);
  return out.str();
}

// Flasher i in metres (ls1..ls4 of ff-gct.cc): z = (356.1075 - 305.03801 + depth) cm.
Vec3 flasher_position(int i, double xy_cm, double depth_cm) {
  const double sx = i < 2 ? 1.0 : -1.0;
  const double sy = i % 2 == 0 ? 1.0 : -1.0;
  return {0.01 * sx * xy_cm, 0.01 * sy * xy_cm, 0.01 * (356.1075 - 305.03801 + depth_cm)};
}

void default_geometry_and_counts() {
  ff_gct::Options options;
  options.photons = {1e5};
  const auto s = read_stream(stream_bytes(options));

  // Four flashers fire in one event: 4 x 1e5 photons (no target inside the fiducial sphere).
  CHECK(count_type(s, gstrickd::RecordType::event) == 1);
  CHECK(count_type(s, gstrickd::RecordType::telescope) == 1);
  CHECK(s.rows.size() == 4 * 100000);
  for (const Row &r : s.rows) {
    CHECK(r.frame == 2); // telescope frame
    CHECK(!r.has_zem);
    CHECK(r.w == 1.0F);
    CHECK(r.t == 0.0F);
    CHECK(r.wavelength == 400.0F);
    CHECK_NEAR(r.cx * r.cx + r.cy * r.cy + r.cz * r.cz, 1.0, 1e-5);
  }
  // Positions of the first flasher (rows 0..99999): (+14.62, +14.62, 53.53) cm, exact.
  const Vec3 p0 = flasher_position(0, 14.62, 2.46);
  for (std::size_t i = 0; i < 100000; ++i) {
    CHECK_NEAR(s.rows[i].x, p0.x, 1e-5);
    CHECK_NEAR(s.rows[i].y, p0.y, 1e-5);
    CHECK_NEAR(s.rows[i].z, p0.z, 1e-5);
  }
  // Isotropic emission: the mean of cz is zero (variance 1/3).
  double sum_cz = 0.0;
  for (std::size_t i = 0; i < 100000; ++i) {
    sum_cz += s.rows[i].cz;
  }
  CHECK_NEAR(sum_cz / 100000.0, 0.0, 4.0 * std::sqrt(1.0 / 3.0 / 100000.0));
  // The last flasher sits at (-14.62, -14.62).
  const Vec3 p3 = flasher_position(3, 14.62, 2.46);
  CHECK_NEAR(s.rows[3 * 100000].x, p3.x, 1e-5);
  CHECK_NEAR(s.rows[3 * 100000].y, p3.y, 1e-5);
}

void delays_and_fire_pattern() {
  // Only flasher 3 fires, with delay delays[2].
  ff_gct::Options options;
  options.photons = {1e4};
  options.fire = "3";
  options.delays = {5.0, 6.0, 7.0, 8.0, 9.0};
  const auto s = read_stream(stream_bytes(options));
  CHECK(s.rows.size() == 10000);
  const Vec3 p = flasher_position(2, 14.62, 2.46);
  for (const Row &r : s.rows) {
    CHECK_NEAR(r.t, 7.0, 1e-5);
    CHECK_NEAR(r.x, p.x, 1e-5);
    CHECK_NEAR(r.y, p.y, 1e-5);
  }

  // Demo: four individual flashes with delays 0..3, then all four at delays[4].
  ff_gct::Options demo;
  demo.photons = {1e3};
  demo.fire = "demo";
  demo.delays = {0.0, 0.0, 0.0, 0.0, 0.0}; // zero delays select the demo values
  const auto d = read_stream(stream_bytes(demo));
  CHECK(d.rows.size() == 8 * 1000);
  CHECK_NEAR(d.rows[0].t, 0.0, 1e-5);
  CHECK_NEAR(d.rows[1000].t, 10.0, 1e-5);
  CHECK_NEAR(d.rows[2000].t, 20.0, 1e-5);
  CHECK_NEAR(d.rows[3000].t, 30.0, 1e-5);
  CHECK_NEAR(d.rows[4000].t, 50.0, 1e-5);
  CHECK_NEAR(d.rows[7999].t, 50.0, 1e-5);

  // Alternate: event k fires flasher k mod 4 only.
  ff_gct::Options alt;
  alt.photons = {1e3};
  alt.events = 4;
  alt.fire = "alternate";
  const auto a = read_stream(stream_bytes(alt));
  CHECK(count_type(a, gstrickd::RecordType::event) == 4);
  CHECK(a.rows.size() == 4 * 1000);
  for (std::size_t e = 0; e < 4; ++e) {
    const Vec3 q = flasher_position(static_cast<int>(e), 14.62, 2.46);
    CHECK_NEAR(a.rows[e * 1000].x, q.x, 1e-5);
    CHECK_NEAR(a.rows[e * 1000].y, q.y, 1e-5);
    CHECK(a.rows[e * 1000].event == static_cast<std::int64_t>(e + 1));
  }
}

} // namespace

int main() {
  return run_checks([] {
    default_geometry_and_counts();
    delays_and_fire_pattern();
  });
}
