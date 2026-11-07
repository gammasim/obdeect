#include "check.hpp"

#include "dlufd/profile.hpp"
#include "dlufd/refraction.hpp"
#include "dlufd/timing.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#ifndef DLUFD_TEST_DATA_DIR
#error "DLUFD_TEST_DATA_DIR must be defined by CMake"
#endif

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kObservationM = 2180.;

struct Photon {
  double cx, cy, x, y, t;
};

// Photon emitted at height h with zenith angle theta (degrees) in the x-z plane.
Photon start(double zenith_deg) {
  const double s = std::sin(zenith_deg * kPi / 180.);
  return {s, 0., 0., 0., 0.};
}

} // namespace

int main() {
  return run_checks([] {
    const auto profile = dlufd::AtmosphereProfile::from_file(
        std::string(DLUFD_TEST_DATA_DIR) + "/atmprof_ecmwf_north_winter_fixed.dat");
    const dlufd::RefractionTables tables(profile, kObservationM);
    CHECK(tables.etadsn() > 0.22 && tables.etadsn() < 0.24);

    // An observation level outside the table is rejected.
    bool rejected = false;
    try {
      const dlufd::RefractionTables bad(profile, 130000.);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    CHECK(rejected);

    // Vertical photon: no displacement or direction change; the time delay equals the
    // straight-path air delay of the profile (both use the same table).
    for (double h : {3000., 10000., 20000.}) {
      Photon p = start(0.);
      const bool bent = dlufd::bend(tables, h, p.cx, p.cy, p.x, p.y, p.t);
      CHECK(bent);
      CHECK(p.cx == 0. && p.cy == 0.);
      CHECK(p.x == 0. && p.y == 0.);
      const double vacuum_ns = (h - kObservationM) / dlufd::speed_of_light_m_per_ns;
      const double air_ns = dlufd::straight_path_air_delay_ns(profile, kObservationM, h, 0.);
      CHECK(p.t > 0.);
      CHECK_NEAR(p.t - vacuum_ns, air_ns, 1e-6);
      CHECK(air_ns > 0.);
    }

    // Displacement grows with zenith angle, is directed against the photon's horizontal
    // direction, and is of order metres at 20 km; the direction stays a unit vector.
    double previous_x = 0.;
    for (double zenith : {20., 40., 60.}) {
      Photon p = start(zenith);
      const double cx_in = p.cx;
      const bool bent = dlufd::bend(tables, 20000., p.cx, p.cy, p.x, p.y, p.t);
      CHECK(bent);
      CHECK(p.x < 0.);
      CHECK(std::fabs(p.x) > std::fabs(previous_x));
      previous_x = p.x;
      CHECK(std::fabs(p.x) < 50.);
      // Refraction bends towards the vertical: the horizontal direction cosine shrinks.
      CHECK(p.cx > 0. && p.cx < cx_in);
      const double cz = std::sqrt(1. - p.cx * p.cx - p.cy * p.cy);
      CHECK_NEAR(p.cx * p.cx + p.cy * p.cy + cz * cz, 1., 1e-12);
    }

    // A downward photon from above the observation level is never totally reflected (n-1
    // decreases with height), so a 75 degree photon is bent, not left unchanged.
    Photon low = start(75.);
    CHECK(dlufd::bend(tables, 20000., low.cx, low.cy, low.x, low.y, low.t));
    CHECK(low.x < 0. && low.cx > 0. && low.cx < std::sin(75. * kPi / 180.));

    // A direction beyond the unit circle (upward or horizontal) is not handled.
    Photon upward{1.1, 0., 0., 0., 0.};
    CHECK(!dlufd::bend(tables, 20000., upward.cx, upward.cy, upward.x, upward.y, upward.t));

    // The float instantiation agrees with double to single precision.
    Photon pd = start(40.);
    Photon pf = start(40.);
    float fcx = static_cast<float>(pf.cx), fcy = 0.f, fx = 0.f, fy = 0.f, ft = 0.f;
    CHECK(dlufd::bend(tables, 15000., pd.cx, pd.cy, pd.x, pd.y, pd.t));
    CHECK(dlufd::bend(tables, 15000.f, fcx, fcy, fx, fy, ft));
    CHECK_NEAR(static_cast<double>(fcx), pd.cx, 1e-6);
    CHECK_NEAR(static_cast<double>(fx), pd.x, 1e-3);
    CHECK_NEAR(static_cast<double>(ft), pd.t, 1e-2);
  });
}
