#include "check.hpp"

#include "dlamb/programs/ls_geo_test.hpp"

#include <cmath>
#include <sstream>
#include <string>

using namespace dlamb;
using namespace dlamb::ls_geo_test;

namespace {

void distances() {
  const Vec3 o{0, 0, 0}, px{10, 0, 0}, py{0, 10, 0}, pz{0, 0, 10};
  CHECK_NEAR(norm(px - o), 10.0, 1e-12);
  CHECK_NEAR(distance(line(px, {0, 0, 1}), line(o, {0, 1, 0})), 10.0, 1e-12);
  CHECK_NEAR(distance(line(py, {1, 0, 0}), line(o, {0, 0, 1})), 10.0, 1e-12);
  CHECK_NEAR(distance(line(pz, {0, 1, 0}), line(o, {1, 0, 0})), 10.0, 1e-12);
  // Parallel lines: point-to-line distance.
  CHECK_NEAR(distance(line(px, {0, 0, 1}), line(o, {0, 0, 5})), 10.0, 1e-12);
  const Vec3 p345{3, 4, 5};
  CHECK_NEAR(distance(p345, line(o, {1, 0, 0})), std::sqrt(41.0), 1e-12);
  CHECK_NEAR(distance(p345, plane(o, {0, 0, 1})), 5.0, 1e-12);
  CHECK_NEAR(distance(o, plane(p345, {1, 2, 0})), 11.0 / std::sqrt(5.0), 1e-12);
  CHECK_NEAR(distance(p345, plane(p345, {1, 2, 0})), 0.0, 1e-12);
}

void rotations_and_intersections() {
  const Vec3 p345{3, 4, 5};
  const double half_pi = pi / 2.0;
  const Vec3 rx = rotate_x(p345, half_pi), ry = rotate_y(p345, half_pi),
             rz = rotate_z(p345, half_pi);
  CHECK_NEAR(rx.x, 3.0, 1e-12);
  CHECK_NEAR(rx.y, -5.0, 1e-12);
  CHECK_NEAR(rx.z, 4.0, 1e-12);
  CHECK_NEAR(ry.x, 5.0, 1e-12);
  CHECK_NEAR(ry.y, 4.0, 1e-12);
  CHECK_NEAR(ry.z, -3.0, 1e-12);
  CHECK_NEAR(rz.x, -4.0, 1e-12);
  CHECK_NEAR(rz.y, 3.0, 1e-12);
  CHECK_NEAR(rz.z, 5.0, 1e-12);

  const Plane a = plane(p345, {1, 2, 0});
  const auto on_x = intersection(line({0, 0, 0}, {1, 0, 0}), a);
  CHECK(on_x.has_value());
  CHECK_NEAR(on_x->x, 11.0, 1e-9);
  const auto on_y = intersection(line({0, 0, 0}, {0, 1, 0}), a);
  CHECK(on_y.has_value());
  CHECK_NEAR(on_y->y, 5.5, 1e-9);
  // The z axis is parallel to the plane (normal has no z component).
  CHECK(!intersection(line({0, 0, 0}, {0, 0, 1}), a).has_value());
}

void report() {
  std::ostringstream os;
  write_report(os);
  const std::string text = os.str();
  CHECK(text.find("Distance O - Px = 10 m") != std::string::npos);
  CHECK(text.find("A has no intersection with z axis") != std::string::npos);
}

} // namespace

int main() {
  return run_checks([] {
    distances();
    rotations_and_intersections();
    report();
  });
}
