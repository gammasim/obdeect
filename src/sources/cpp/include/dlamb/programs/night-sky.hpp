#pragma once

// Night-sky background photons on one telescope sphere, sky frame (docs/spec/tool-night-sky.md).

#include "dlamb/sky_sources.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>

namespace dlamb::programs {

struct NsbOptions {
  double telescope_radius_m{15.0};
  double viewcone_deg{10.0}; // half opening angle
  double altitude_m{2150.0};
  std::int64_t atmosphere{-1};
  double photons{1e6}; // per event
  double bunch_size{2.0};
  std::uint32_t events{1};
  std::int64_t run_number{1};
  Spectrum spectrum{Fixed{400.0}};
  double pulse_width_ns{0.0};
  double pulse_delay_ns{0.0};
  std::uint64_t seed{1};
};

inline SkyPlan<BackgroundEmission> plan_night_sky(const NsbOptions &o) {
  if (!std::isfinite(o.telescope_radius_m) || !(o.telescope_radius_m > 0.0))
    throw std::invalid_argument("telescope radius must be finite and > 0");
  if (!std::isfinite(o.viewcone_deg) || !(o.viewcone_deg > 0.0 && o.viewcone_deg <= 180.0))
    throw std::invalid_argument("view cone must lie in (0, 180] deg");
  if (!std::isfinite(o.altitude_m))
    throw std::invalid_argument("observation level must be finite");
  if (o.atmosphere < -1)
    throw std::invalid_argument("atmosphere: profile number (or -1 for none)");
  if (!std::isfinite(o.photons) || o.photons < 0.0 || o.photons > max_photon_count)
    throw std::invalid_argument("photons must be finite, >= 0 and at most 2^53");
  if (!std::isfinite(o.bunch_size) || !(o.bunch_size > 0.0))
    throw std::invalid_argument("bunch size must be finite and > 0");
  if (!std::isfinite(o.pulse_width_ns) || o.pulse_width_ns < 0.0)
    throw std::invalid_argument("pulse width must be finite and >= 0");
  if (!std::isfinite(o.pulse_delay_ns))
    throw std::invalid_argument("pulse delay must be finite");
  validate_spectrum(o.spectrum);

  const double bunch = std::max(1.0, o.bunch_size);
  SkyPlan<BackgroundEmission> plan;
  plan.run.run_number = o.run_number > 0 ? o.run_number : 1;
  plan.run.observation_level_m = o.altitude_m;
  plan.run.atmosphere_id = o.atmosphere;
  const auto [lo, hi] = spectrum_range(o.spectrum);
  plan.run.wavelength_min_nm = lo;
  plan.run.wavelength_max_nm = hi;
  plan.run.telescopes.push_back({1, {0.0, 0.0, 0.0}, o.telescope_radius_m});
  plan.run.extra["viewcone_deg"] = o.viewcone_deg;
  plan.run.extra["background"] = std::string("night sky");
  plan.events_per_case = o.events;

  BackgroundEmission e;
  e.source.centre = {0.0, 0.0, 0.0};
  e.source.radius_m = o.telescope_radius_m;
  e.source.pointing = {0.0, 0.0, -1.0};
  e.source.viewcone_rad = o.viewcone_deg * pi / 180.0;
  e.source.spectrum = o.spectrum;
  e.source.pulse_width_ns = o.pulse_width_ns;
  e.source.pulse_delay_ns = o.pulse_delay_ns;
  e.seed = o.seed;
  e.bunch = bunch;
  e.mean_bunches = o.photons / bunch;
  plan.cases.push_back({e});
  return plan;
}

inline void write_night_sky(std::ostream &out, const NsbOptions &o) {
  write_plan(out, "dlamb-night-sky", plan_night_sky(o));
}

} // namespace dlamb::programs
