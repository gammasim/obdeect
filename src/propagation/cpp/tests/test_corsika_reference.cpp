// Agreement with a reference implementation of the profile interpolation and ray bending, built
// only with -DDLUFD_CORSIKA_REFERENCE=<reference source directory>. Same profile file, same
// observation level.
#include "check.hpp"

#include "dlufd/profile.hpp"
#include "dlufd/refraction.hpp"

#include <cmath>
#include <string>

#ifndef DLUFD_TEST_DATA_DIR
#error "DLUFD_TEST_DATA_DIR must be defined by CMake"
#endif

extern "C" {
void atmnam_(const char *aname, double *obslev);
double rhofx_(double *height);
double thickx_(double *height);
double refim1x_(double *height);
void raybnd_(double *zem, double *u, double *v, double *w, double *dx, double *dy, double *dt);

// Built-in-atmosphere functions the reference code references; unused with a profile file.
double rhof_(double *) { return 0.; }
double heigh_(double *) { return 0.; }
double thick_(double *) { return 0.; }
}

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kObservationM = 2180.;
} // namespace

int main() {
  return run_checks([] {
    const std::string path =
        std::string(DLUFD_TEST_DATA_DIR) + "/atmprof_ecmwf_north_winter_fixed.dat";
    double obs_cm = 100. * kObservationM;
    atmnam_(path.c_str(), &obs_cm);

    const auto profile = dlufd::AtmosphereProfile::from_file(path);
    const dlufd::RefractionTables tables(profile, kObservationM);

    // Lookups in cm: the same linear interpolation on the same table.
    double worst_rho = 0., worst_thick = 0., worst_n1 = 0.;
    for (double h = 0.; h < 119.e5; h += 0.0137e5) {
      double height = h;
      const double rho = rhofx_(&height), thick = thickx_(&height), n1 = refim1x_(&height);
      worst_rho = std::fmax(worst_rho, std::fabs(profile.density_cm(h) - rho) / rho);
      if (thick > 0.)
        worst_thick = std::fmax(worst_thick, std::fabs(profile.thickness_cm(h) - thick) / thick);
      if (n1 > 0.)
        worst_n1 = std::fmax(worst_n1, std::fabs(profile.index_minus_one_cm(h) - n1) / n1);
    }
    CHECK(worst_rho < 1e-12);
    CHECK(worst_thick < 1e-12);
    CHECK(worst_n1 < 1e-12);

    // Ray bending: displacement [cm -> m], direction and time, for several emission heights.
    for (double height_m : {2500., 5000., 10000., 20000., 40000.}) {
      for (double zenith_deg : {5., 20., 40., 60.}) {
        const double u = std::sin(zenith_deg * kPi / 180.);
        const double v = 0.3 * u;
        const double w = std::sqrt(1. - u * u - v * v);
        double zem = 100. * height_m, cu = u, cv = v, cw = w, dx = 0., dy = 0., dt = 0.;
        raybnd_(&zem, &cu, &cv, &cw, &dx, &dy, &dt);

        double cx = u, cy = v, x = 0., y = 0., t = 0.;
        const bool bent = dlufd::bend(tables, height_m, cx, cy, x, y, t);
        CHECK(bent);
        CHECK_NEAR(cx, cu, 1e-12);
        CHECK_NEAR(cy, cv, 1e-12);
        CHECK_NEAR(x, dx / 100., 1e-9);
        CHECK_NEAR(y, dy / 100., 1e-9);
        CHECK_NEAR(t, dt, 1e-6);
      }
    }

    // Vertical emission: the reference adds the travel time only.
    double zem = 1000000., cu = 0., cv = 0., cw = 1., dx = 0., dy = 0., dt = 0.;
    raybnd_(&zem, &cu, &cv, &cw, &dx, &dy, &dt);
    double cx = 0., cy = 0., x = 0., y = 0., t = 0.;
    CHECK(dlufd::bend(tables, 10000., cx, cy, x, y, t));
    CHECK_NEAR(t, dt, 1e-6);
  });
}
