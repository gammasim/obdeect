#pragma once

// Port of sim_telarray/LightEmission/xyls.cc: a parallel (vertical, downward) light beam
// scanned along two straight lines over the camera area. Test program of the original,
// kept for comparison. Output: sky-frame photons, one emission per scan position.

#include "dlamb/sky_sources.hpp"
#include "dlamb/source.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace dlamb::programs {

// Options of xyls.cc; units cm unless noted. Defaults of x, y and height are derived from
// radius 123 cm in the original (as its declarations do, not from a later -r).
struct XylsOptions {
  double radius_cm{123.0};            // -r, camera fiducial radius (also scan length)
  double x_cm{0.3 * 123.0};           // -x, start of the scan (default 0.3 r of 123 cm)
  double y_cm{0.2 * 123.0};           // -y, start of the scan (default 0.2 r of 123 cm)
  double height_cm{0.05 * 123.0};     // -h, height of the source (default 0.05 r of 123 cm)
  double photons_per_position{100.0}; // nphot_per_pos (fixed in xyls.cc)
  double bunch_size{1.0};             // bunch_size (fixed in xyls.cc)
  double altitude_m{2150.0};          // alt_m (fixed in xyls.cc)
  std::int64_t atmosphere{1};         // Run(1, ..., 1, ...): recorded only
  std::uint64_t seed{1};              // random seed (new option)
};

[[nodiscard]] inline SkyPlan<SourceEmission> plan_xyls(const XylsOptions &o) {
  SkyPlan<SourceEmission> plan;
  plan.run.run_number = 1;
  plan.run.observation_level_m = o.altitude_m;
  plan.run.atmosphere_id = o.atmosphere;
  plan.run.wavelength_min_nm = 400.0;
  plan.run.wavelength_max_nm = 400.0;
  plan.run.telescopes.push_back({1, {0.0, 0.0, 0.0}, o.radius_cm / 100.0});
  plan.events_per_case = 1;

  const double bunch = std::max(1.0, o.bunch_size);
  const double r = o.radius_cm;
  if (!(o.photons_per_position >= 0.0) || o.photons_per_position > max_photon_count)
    throw std::invalid_argument("xyls: photons per position must be finite, >= 0 and at most 2^53");
  std::vector<SourceEmission> list;
  const auto add = [&](double x_cm, double y_cm) {
    Source s;
    s.frame = gstrickd::Frame::sky;
    s.geometry = Emitter{
        {x_cm / 100.0, y_cm / 100.0, o.height_cm / 100.0}, {0.0, 0.0, -1.0}, FixedAngle{0.0}};
    s.spectrum = Fixed{400.0};
    s.pulse = GaussianShape{0.0, 0.0};
    s.photons = o.photons_per_position / bunch;
    list.push_back({s, o.seed + list.size(), bunch, 0.0});
  };
  // Same loops as xyls.cc: a line along x at y = y0, then a line along y at x = x0.
  for (double x = o.x_cm; x < o.x_cm + 0.4 * r; x += 0.001 * r)
    add(x, o.y_cm);
  for (double y = o.y_cm; y < o.y_cm + 0.2 * r; y += 0.001 * r)
    add(o.x_cm, y);
  plan.cases.push_back(std::move(list));
  return plan;
}

inline void write_xyls(std::ostream &out, const XylsOptions &o) {
  write_plan(out, "dlamb-xyls", plan_xyls(o));
}

} // namespace dlamb::programs
