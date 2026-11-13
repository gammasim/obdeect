#pragma once

// Point light source above the observation level, sky frame (docs/spec/tool-point-source.md).

#include "obdeect/sources/sky_sources.hpp"
#include "obdeect/stream/telpos.hpp"

#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace obdeect::sources::programs {

using obdeect::stream::cm_to_m;
using obdeect::stream::TelposCm;

struct XyzlsOptions {
  std::int64_t run_number{1};
  std::uint32_t events{1}; // per case
  double telcam_radius_cm{1500.0};
  double x_cm{0.0}, y_cm{0.0};
  double z_cm{1000e2};
  std::vector<double> photons{1e7}; // into 4 pi per event; one case per value
  double bunch_size{1.0};
  double altitude_m{2150.0};
  std::int64_t atmosphere{1};
  Vec3 direction{0.0, 0.0, -1.0};
  Spectrum spectrum{Fixed{400.0}};
  Pulse pulse{GaussianShape{0.0, 0.0}};
  Angular angular{Isotropic{}};
  std::uint64_t seed{1};
  std::vector<TelposCm> telescopes; // empty: one telescope of radius telcam_radius_cm
};

inline SkyPlan<SourceEmission> plan_point_source(const XyzlsOptions &o) {
  if (!within_limits(o.x_cm) || !within_limits(o.y_cm) || !within_limits(o.z_cm))
    throw std::invalid_argument("source position out of range");
  if (!std::isfinite(o.altitude_m))
    throw std::invalid_argument("observation level must be finite");
  if (!std::isfinite(o.bunch_size))
    throw std::invalid_argument("bunch size must be finite");
  if (!is_finite(o.direction) || norm(o.direction) == 0.0)
    throw std::invalid_argument("direction vector must not be zero");
  SkyPlan<SourceEmission> plan;
  plan.run.run_number = o.run_number > 0 ? o.run_number : 1;
  plan.run.observation_level_m = o.altitude_m;
  plan.run.atmosphere_id = o.atmosphere;
  const auto [lo, hi] = spectrum_range(o.spectrum);
  plan.run.wavelength_min_nm = lo;
  plan.run.wavelength_max_nm = hi;
  if (o.telescopes.empty()) {
    if (!(o.telcam_radius_cm > 0.0) || !std::isfinite(o.telcam_radius_cm))
      throw std::invalid_argument("telescope radius must be finite and positive");
    plan.run.telescopes.push_back({1, {0.0, 0.0, 0.0}, o.telcam_radius_cm / 100.0});
  } else {
    std::uint32_t id = 1;
    for (const auto &t : o.telescopes)
      plan.run.telescopes.push_back(
          {id++, {cm_to_m(t.x), cm_to_m(t.y), cm_to_m(t.z)}, cm_to_m(t.r)});
  }
  plan.events_per_case = o.events;

  const double bunch = std::max(1.0, o.bunch_size);
  const Vec3 position{o.x_cm / 100.0, o.y_cm / 100.0, o.z_cm / 100.0};
  const double rmax = fiducial_super_radius(plan.run.telescopes);
  Target target;
  if (norm(position) > 1.02 * rmax) {
    if (std::holds_alternative<Isotropic>(o.angular))
      target = SphereTarget{{0.0, 0.0, 0.0}, 1.02 * rmax};
    else if (!std::holds_alternative<FixedAngle>(o.angular))
      throw std::invalid_argument("a non-isotropic profile beyond 1.02 x the fiducial radius is "
                                  "not supported (use isotropic or Cone/Parallel)");
  }
  for (std::size_t c = 0; c < o.photons.size(); ++c) {
    const double n = o.photons[c];
    if (!std::isfinite(n) || n < 0.0 || n > max_photon_count)
      throw std::invalid_argument("photons must be finite, >= 0 and at most 2^53");
    SourceEmission e;
    e.source.frame = obdeect::stream::Frame::sky;
    e.source.geometry = Emitter{position, normalised(o.direction), o.angular};
    e.source.spectrum = o.spectrum;
    e.source.pulse = o.pulse;
    e.source.photons = n / bunch;
    e.source.target = target;
    validate(e.source);
    e.seed = o.seed + c;
    e.bunch = bunch;
    plan.cases.push_back({e});
  }
  return plan;
}

inline void write_point_source(std::ostream &out, const XyzlsOptions &o) {
  write_plan(out, "obdeect sources point-source", plan_point_source(o));
}

} // namespace obdeect::sources::programs
