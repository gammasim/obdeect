#pragma once

// Parallel downward beam scanned along two lines over the camera area, sky frame
// (docs/spec/tool-beam-scan.md).

#include "dlamb/sky_sources.hpp"

#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>

namespace dlamb::programs {

struct XylsOptions {
  double radius_cm{123.0};
  double x_cm{0.3 * 123.0};
  double y_cm{0.2 * 123.0};
  double height_cm{0.05 * 123.0};
  double photons_per_position{100.0};
  double bunch_size{1.0};
  std::uint64_t seed{1};
};

inline SkyPlan<SourceEmission> plan_beam_scan(const XylsOptions &o) {
  if (!std::isfinite(o.radius_cm) || !(o.radius_cm > 0.0))
    throw std::invalid_argument("radius must be finite and positive");
  if (!within_limits(o.x_cm) || !within_limits(o.y_cm) || !within_limits(o.height_cm))
    throw std::invalid_argument("scan start or height out of range");
  if (!std::isfinite(o.photons_per_position) || o.photons_per_position < 0.0 ||
      o.photons_per_position > max_photon_count)
    throw std::invalid_argument("photons per position must be finite, >= 0 and at most 2^53");
  if (!std::isfinite(o.bunch_size))
    throw std::invalid_argument("bunch size must be finite");
  const double bunch = std::max(1.0, o.bunch_size);
  SkyPlan<SourceEmission> plan;
  plan.run.run_number = 1;
  plan.run.observation_level_m = 2150.0;
  plan.run.atmosphere_id = 1;
  plan.run.wavelength_min_nm = 400.0;
  plan.run.wavelength_max_nm = 400.0;
  plan.run.telescopes.push_back({1, {0.0, 0.0, 0.0}, o.radius_cm / 100.0});
  plan.events_per_case = 1;

  std::vector<SourceEmission> scan;
  const auto add = [&](double x_cm, double y_cm) {
    SourceEmission e;
    e.source.frame = gstrickd::Frame::sky;
    e.source.geometry = Emitter{
        {x_cm / 100.0, y_cm / 100.0, o.height_cm / 100.0}, {0.0, 0.0, -1.0}, FixedAngle{0.0}};
    e.source.spectrum = Fixed{400.0};
    e.source.pulse = GaussianShape{0.0, 0.0};
    e.source.photons = o.photons_per_position / bunch;
    validate(e.source);
    e.seed = o.seed + scan.size();
    e.bunch = bunch;
    scan.push_back(e);
  };
  for (double x = o.x_cm; x < o.x_cm + 0.4 * o.radius_cm; x += 0.001 * o.radius_cm)
    add(x, o.y_cm);
  for (double y = o.y_cm; y < o.y_cm + 0.2 * o.radius_cm; y += 0.001 * o.radius_cm)
    add(o.x_cm, y);
  plan.cases.push_back(std::move(scan));
  return plan;
}

inline void write_beam_scan(std::ostream &out, const XylsOptions &o) {
  write_plan(out, "dlamb-beam-scan", plan_beam_scan(o));
}

} // namespace dlamb::programs
