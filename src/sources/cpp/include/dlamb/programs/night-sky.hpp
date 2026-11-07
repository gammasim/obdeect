#pragma once

// Port of sim_telarray/LightEmission/nsbls.cc: night-sky background photons on one
// telescope sphere, in a viewing cone around the zenith (downward directions).

#include "dlamb/background.hpp"
#include "dlamb/sky_sources.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace dlamb::programs {

// Options of nsbls.cc (units as in the original: m, ns, deg, nm).
struct NsbOptions {
  double telescope_radius_m{15.0}; // -r (--radius, --telescope-radius)
  double viewcone_deg{10.0};       // -v (--viewcone): half opening angle
  double altitude_m{2150.0};       // --altitude
  std::int64_t atmosphere{-1};     // --atmosphere: recorded only (original: none)
  double photons{1e6};             // --photons, photons per event (bunch counting)
  double bunch_size{2.0};          // -b (values below 1 are used as 1)
  std::uint32_t events{1};         // --events
  std::int64_t run_number{1};      // --run
  // -s: the original defaults to Benn_LaPalma_sky_converted.lis (sim_telarray/cfg/common).
  // The default here is 400 nm; the tool reads the table when given (or by default).
  Spectrum spectrum{Fixed{400.0}};
  double pulse_width_ns{0.0}; // -w (--pulse-width): see sample_background
  double pulse_delay_ns{0.0}; // -d (--pulse-delay)
  std::uint64_t seed{1};      // random seed (new option)
};

struct NsbEmission {
  BackgroundSource source{};
  std::uint64_t seed{};
  double bunch{1.0};
  double mean_bunches{};

  [[nodiscard]] std::uint64_t count(std::uint32_t event) const {
    return background_count(mean_bunches, seed, event, shared_telescope_key);
  }
  void fill(std::uint32_t event, std::uint64_t first, std::size_t rows,
            gstrickd::PhotonChunk &chunk) const {
    emit_background(source, seed, event, shared_telescope_key, first, rows, bunch, chunk);
  }
};

[[nodiscard]] inline SkyPlan<NsbEmission> plan_night_sky(const NsbOptions &o) {
  const double bunch = std::max(1.0, o.bunch_size);
  const auto [wl_min, wl_max] = spectrum_range(o.spectrum);

  SkyPlan<NsbEmission> plan;
  plan.run.run_number = o.run_number > 0 ? o.run_number : 1;
  plan.run.observation_level_m = o.altitude_m;
  plan.run.atmosphere_id = o.atmosphere;
  plan.run.wavelength_min_nm = wl_min;
  plan.run.wavelength_max_nm = wl_max;
  plan.run.telescopes.push_back({1, {0.0, 0.0, 0.0}, o.telescope_radius_m});
  plan.run.extra["viewcone_deg"] = o.viewcone_deg;
  plan.run.extra["background"] = std::string("night sky");
  plan.events_per_case = o.events;

  BackgroundSource source;
  source.centre = {0.0, 0.0, 0.0};
  source.radius_m = o.telescope_radius_m;
  source.pointing = {0.0, 0.0, -1.0};
  source.viewcone_rad = o.viewcone_deg * pi / 180.0;
  source.spectrum = o.spectrum;
  source.pulse_width_ns = o.pulse_width_ns;
  source.pulse_delay_ns = o.pulse_delay_ns;
  if (!(o.photons >= 0.0) || o.photons > max_photon_count)
    throw std::invalid_argument("night-sky: photons must be finite, >= 0 and at most 2^53");
  plan.cases.push_back({NsbEmission{source, o.seed, bunch, o.photons / bunch}});
  return plan;
}

inline void write_night_sky(std::ostream &out, const NsbOptions &o) {
  write_plan(out, "dlamb-night-sky", plan_night_sky(o));
}

} // namespace dlamb::programs
