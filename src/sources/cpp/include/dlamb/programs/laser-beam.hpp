#pragma once

// Primary photons of a laser beam inside the field-of-view cone of a telescope, sky frame
// (docs/spec/tool-laser-beam.md). Scattering is not generated.

#include "dlamb/sky_sources.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>

namespace dlamb::programs {

struct LsBeamOptions {
  std::uint32_t events{1};
  double step_m{0.1};
  double bunch_size{3.0};
  double primary_photons_per_m{1e4};
  Spectrum spectrum{Fixed{355.0}};
  Pulse pulse{Fixed{0.0}};
  std::string angular{"Rayleigh"}; // stored verbatim in the run record
  std::int64_t run_number{1};
  double altitude_m{2150.0};
  std::int64_t atmosphere{1};
  Vec3 laser_position_m{1000.0, 0.0, 0.0};
  double laser_theta_deg{60.0};
  double laser_phi_deg{180.0}; // mathematical azimuth
  double telescope_theta_deg{30.0};
  double telescope_phi_deg{0.0}; // mathematical azimuth
  double telescope_radius_m{12.5};
  double fov_deg{15.0};
  std::uint64_t seed{1};
};

// Where the beam enters the field-of-view cone (relative to the telescope), its unit direction
// and the path length inside the cone.
struct BeamGeometry {
  Vec3 start_m;
  Vec3 direction;
  double path_m{};
};

inline BeamGeometry laser_beam_geometry(const LsBeamOptions &o) {
  constexpr double deg = pi / 180.0;
  if (!std::isfinite(o.fov_deg) || !(o.fov_deg > 0.0 && o.fov_deg < 90.0))
    throw std::invalid_argument("field of view must lie in (0, 90) deg");
  if (!is_finite(o.laser_position_m) || !std::isfinite(o.laser_theta_deg) ||
      !std::isfinite(o.laser_phi_deg) || !std::isfinite(o.telescope_theta_deg) ||
      !std::isfinite(o.telescope_phi_deg))
    throw std::invalid_argument("laser and telescope pointing must be finite");
  const Vec3 dir = direction_from_altitude((90.0 - o.laser_theta_deg) * deg, o.laser_phi_deg * deg);
  const double phi_t = o.telescope_phi_deg * deg, theta_t = o.telescope_theta_deg * deg;
  const Vec3 p = rotate_y(rotate_z(o.laser_position_m, -phi_t), -theta_t);
  const Vec3 q = rotate_y(rotate_z(dir, -phi_t), -theta_t);
  const double tan_half = std::tan(0.5 * o.fov_deg * deg);
  const double t2 = tan_half * tan_half;
  const double a = q.x * q.x + q.y * q.y - q.z * q.z * t2;
  const double b = 2.0 * (q.x * p.x + q.y * p.y - q.z * p.z * t2);
  const double c = p.x * p.x + p.y * p.y - p.z * p.z * t2;
  const char *miss = "laser beam does not cross the field of view";
  double r1, r2;
  if (a == 0.0) {
    if (b == 0.0)
      throw std::invalid_argument(miss);
    r1 = r2 = -c / b;
  } else {
    const double disc = b * b - 4.0 * a * c;
    if (disc < 0.0)
      throw std::invalid_argument(miss);
    // q = -(b + sign(b) sqrt(disc)) / 2 avoids cancellation between the two roots.
    const double s = std::sqrt(disc);
    const double qq = -0.5 * (b + (b >= 0.0 ? s : -s));
    if (qq == 0.0) {
      r1 = r2 = 0.0;
    } else {
      r1 = qq / a;
      r2 = c / qq;
    }
  }
  if (r1 > r2)
    std::swap(r1, r2);
  return {o.laser_position_m + dir * r1, dir, r2 - r1};
}

inline SkyPlan<SourceEmission> plan_laser_beam(const LsBeamOptions &o) {
  const BeamGeometry g = laser_beam_geometry(o);
  if (!std::isfinite(o.step_m) || !(o.step_m > 0.0))
    throw std::invalid_argument("step length must be finite and > 0");
  if (!std::isfinite(o.primary_photons_per_m) || o.primary_photons_per_m < 0.0)
    throw std::invalid_argument("primary photons per metre must be finite and >= 0");
  if (!std::isfinite(o.bunch_size) || !std::isfinite(o.altitude_m))
    throw std::invalid_argument("bunch size and altitude must be finite");
  if (!(g.path_m / o.step_m < 1e8))
    throw std::invalid_argument("too many beam steps");
  if (!(o.telescope_radius_m > 0.0) || !std::isfinite(o.telescope_radius_m))
    throw std::invalid_argument("telescope radius must be finite and > 0");
  const double bunch = std::max(1.0, o.bunch_size);

  SkyPlan<SourceEmission> plan;
  plan.run.run_number = o.run_number > 0 ? o.run_number : 1;
  plan.run.observation_level_m = o.altitude_m;
  plan.run.atmosphere_id = o.atmosphere;
  const auto [lo, hi] = spectrum_range(o.spectrum);
  plan.run.wavelength_min_nm = lo;
  plan.run.wavelength_max_nm = hi;
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
  double s = 0.0;
  for (std::uint64_t k = 0; s < g.path_m; ++k, s += o.step_m) {
    SourceEmission e;
    e.source.frame = gstrickd::Frame::sky;
    e.source.geometry = Emitter{g.start_m + g.direction * s, g.direction, FixedAngle{0.0}};
    e.source.spectrum = o.spectrum;
    e.source.pulse = o.pulse;
    e.source.photons = o.primary_photons_per_m * o.step_m / bunch;
    validate(e.source);
    e.seed = o.seed + k;
    e.bunch = bunch;
    e.time_offset_ns = s / speed_of_light_m_per_ns;
    steps.push_back(e);
  }
  plan.cases.push_back(std::move(steps));
  return plan;
}

inline void write_laser_beam(std::ostream &out, const LsBeamOptions &o) {
  write_plan(out, "dlamb-laser-beam", plan_laser_beam(o));
}

} // namespace dlamb::programs
