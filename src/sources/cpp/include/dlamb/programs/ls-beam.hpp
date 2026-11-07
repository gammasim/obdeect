#pragma once

// Port of sim_telarray/LightEmission/ls-beam.cc, changed for the chain (see
// docs/programs/ls-beam.md): the laser beam's primary photons are emitted along the beam inside the
// field-of-view cone. The light scattered towards the telescope is not generated here: it is
// dlufd's single-scattering effect. The original emits the scattered light itself.

#include "dlamb/sky_sources.hpp"
#include "dlamb/source.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace dlamb::programs {

// Options of ls-beam.cc (units as in the original: m and deg).
struct LsBeamOptions {
  std::uint32_t events{1}; // --events (default 1)
  double step_m{0.1};      // --step (original: steplength_cm = 100 x step)
  double bunch_size{3.0};  // --bunchsize (values below 1 are used as 1)
  // Replaces --brightness (photons into 4 pi per m at ground, scattered). The primary
  // photons per metre are what dlufd scatters: scattered yield per m = this x beta_s(h)
  // (see docs/programs/ls-beam.md). Default is a sampling choice, not a physics value.
  double primary_photons_per_m{1e4};
  Spectrum spectrum{Fixed{355.0}};         // --spectrum (default 355 nm)
  Pulse pulse{Fixed{0.0}};                 // --lightpulse (default Simple:0)
  std::string angular{"Rayleigh"};         // --angular-distribution: recorded for dlufd
  std::int64_t run_number{1};              // --run
  double altitude_m{2150.0};               // --altitude (original help says 2100, code 2150)
  std::int64_t atmosphere{1};              // --atmosphere: recorded only
  Vec3 laser_position_m{1000.0, 0.0, 0.0}; // --laser-position x,y,z (w.r.t. telescope)
  double laser_theta_deg{60.0};            // --laser-theta (zenith angle of the beam)
  double laser_phi_deg{180.0};             // --laser-phi (mathematical azimuth)
  double telescope_theta_deg{30.0};        // --telescope-theta (zenith of the axis)
  double telescope_phi_deg{0.0};           // --telescope-phi (mathematical azimuth)
  double telescope_radius_m{12.5};         // --telescope-radius
  double fov_deg{15.0};                    // --fov (diameter of the cone in the original)
  std::uint64_t seed{1};                   // random seed (new option)
};

// Intersection of the beam with the field-of-view cone (ls-beam.cc, solve_quadratic_equation).
inline bool solve_quadratic(double a, double b, double c, double &x1, double &x2) {
  const double bb = b * b, ac4 = 4.0 * a * c;
  if (a == 0.0) {
    if (b == 0.0)
      return false;
    x1 = x2 = -c / b;
    return true;
  }
  const double d = bb - ac4;
  if (d < 0.0)
    return false;
  const double sq = std::sqrt(d);
  x1 = (-b + sq) / (2.0 * a);
  x2 = (-b - sq) / (2.0 * a);
  if (std::abs(ac4) < 1e-7 * bb) {
    const double r = -b / (2.0 * a), q = c / a, s = q / (r * r);
    if (r > 0.0)
      x2 = -r * s * (0.5 + 0.125 * s);
    else
      x1 = r * s * (0.5 + 0.125 * s);
  }
  if (x1 > x2)
    std::swap(x1, x2);
  return true;
}

struct BeamGeometry {
  Vec3 start_m{};   // first point of the beam inside the cone (sky frame)
  Vec3 direction{}; // unit propagation direction of the beam
  double path_m{};  // length of the beam inside the cone
};

[[nodiscard]] inline BeamGeometry ls_beam_geometry(const LsBeamOptions &o) {
  const double deg = pi / 180.0;
  const Vec3 dir = direction_from_altitude((90.0 - o.laser_theta_deg) * deg, o.laser_phi_deg * deg);
  // Frame of the telescope axis: Rotate_z(-phi), then Rotate_y(-theta), as ls-beam.cc.
  const Vec3 p = rotate_y(rotate_z(o.laser_position_m, -o.telescope_phi_deg * deg),
                          -o.telescope_theta_deg * deg);
  const Vec3 q = rotate_y(rotate_z(dir, -o.telescope_phi_deg * deg), -o.telescope_theta_deg * deg);
  const double thf = std::tan(0.5 * o.fov_deg * deg);
  double t1 = 0.0, t2 = 0.0;
  if (!solve_quadratic(q.x * q.x + q.y * q.y - q.z * q.z * thf * thf,
                       2.0 * (q.x * p.x + q.y * p.y - q.z * p.z * thf * thf),
                       p.x * p.x + p.y * p.y - p.z * p.z * thf * thf, t1, t2))
    throw std::invalid_argument("ls-beam: laser beam does not cross the field of view");
  return {o.laser_position_m + dir * t1, dir, t2 - t1};
}

// The beam is emitted point by point: at every step of `step_m` a point emitter with
// the beam direction (FixedAngle 0) emits primary_photons_per_m x step photons; time is
// the distance along the beam divided by c (dlufd adds the group delay).
[[nodiscard]] inline SkyPlan<SourceEmission> plan_ls_beam(const LsBeamOptions &o) {
  if (!(o.step_m > 0.0) || !std::isfinite(o.step_m))
    throw std::invalid_argument("ls-beam: step must be positive");
  if (!(o.primary_photons_per_m >= 0.0) || !std::isfinite(o.primary_photons_per_m))
    throw std::invalid_argument("ls-beam: primary photons per m must be finite and >= 0");
  const BeamGeometry g = ls_beam_geometry(o);
  const double bunch = std::max(1.0, o.bunch_size);
  const auto [wl_min, wl_max] = spectrum_range(o.spectrum);

  SkyPlan<SourceEmission> plan;
  plan.run.run_number = o.run_number > 0 ? o.run_number : 1;
  plan.run.observation_level_m = o.altitude_m;
  plan.run.atmosphere_id = o.atmosphere;
  plan.run.wavelength_min_nm = wl_min;
  plan.run.wavelength_max_nm = wl_max;
  plan.run.telescopes.push_back({1, {0.0, 0.0, 0.0}, o.telescope_radius_m});
  plan.run.extra["scattering_angular"] = o.angular;
  plan.run.extra["primary_photons_per_m"] = o.primary_photons_per_m;
  plan.run.extra["beam_step_m"] = o.step_m;
  plan.run.extra["beam_path_m"] = g.path_m;
  plan.run.extra["fov_deg"] = o.fov_deg;
  plan.run.extra["telescope_zenith_deg"] = o.telescope_theta_deg;
  plan.run.extra["telescope_azimuth_deg"] = o.telescope_phi_deg;
  plan.events_per_case = o.events;

  std::vector<SourceEmission> steps;
  std::uint64_t k = 0;
  for (double s = 0.0; s < g.path_m; s += o.step_m, ++k) {
    Source src;
    src.frame = gstrickd::Frame::sky;
    src.geometry = Emitter{g.start_m + g.direction * s, g.direction, FixedAngle{0.0}};
    src.spectrum = o.spectrum;
    src.pulse = o.pulse;
    src.photons = o.primary_photons_per_m * o.step_m / bunch;
    steps.push_back({src, o.seed + k, bunch, s / speed_of_light_m_per_ns});
  }
  plan.cases.push_back(std::move(steps));
  return plan;
}

inline void write_ls_beam(std::ostream &out, const LsBeamOptions &o) {
  write_plan(out, "dlamb-ls-beam", plan_ls_beam(o));
}

} // namespace dlamb::programs
