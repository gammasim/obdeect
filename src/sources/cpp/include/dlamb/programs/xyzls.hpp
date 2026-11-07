#pragma once

// Port of sim_telarray/LightEmission/xyzls.cc: a point light source at (x, y, z) above
// the observation level, emitting bunches into 4 pi (or into the fiducial sphere cone).
// Output: sky-frame photons of one telescope block per telescope.

#include "dlamb/sky_sources.hpp"
#include "dlamb/source.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace dlamb::programs {

// Telescope of --telpos-file / --corsika-telpos-file (centre and radius in cm).
struct TelposCm {
  double x{};
  double y{};
  double z{};
  double r{};
};

// Options of xyzls.cc; units as in the original (cm unless noted). Line numbers refer
// to sim_telarray/LightEmission/xyzls.cc.
struct XyzlsOptions {
  std::int64_t run_number{1};           // -R (main, default 1)
  std::uint32_t events{1};              // -N (default 1)
  double telcam_radius_cm{1500.0};      // -r, fiducial radius of the default telescope
  double x_cm{0.0};                     // -x (default 0)
  double y_cm{0.0};                     // -y (default 0)
  double z_cm{1000e2};                  // -z, height above the observation level
  std::vector<double> photons{1e7};     // -n, one case per entry, `events` events each
  double bunch_size{1.0};               // -b (values below 1 are used as 1, as in Emit)
  double altitude_m{2150.0};            // -h, observation level [m a.s.l.]
  std::int64_t atmosphere{1};           // -A, recorded in the run keys (not applied here)
  Vec3 direction{0.0, 0.0, -1.0};       // -d cx,cy,cz (normalised)
  Spectrum spectrum{Fixed{400.0}};      // -s (nm, a table, or Cherenkov:a:b)
  Pulse pulse{GaussianShape{0.0, 0.0}}; // -p (default Gauss:0.0)
  Angular angular{Isotropic{}};         // -a (default isotropic)
  std::uint64_t seed{1};                // random seed (new option)
  std::vector<TelposCm> telescopes{};   // --telpos-file; empty: one telescope at (0,0,0)
};

[[nodiscard]] inline SkyRun xyzls_run(const XyzlsOptions &o) {
  SkyRun run;
  run.run_number = o.run_number > 0 ? o.run_number : 1; // Run::init
  run.observation_level_m = o.altitude_m;
  run.atmosphere_id = o.atmosphere;
  const auto [wl_min, wl_max] = spectrum_range(o.spectrum);
  run.wavelength_min_nm = wl_min;
  run.wavelength_max_nm = wl_max;
  if (o.telescopes.empty()) {
    run.telescopes.push_back({1, {0.0, 0.0, 0.0}, o.telcam_radius_cm / 100.0});
  } else {
    std::uint32_t id = 0;
    for (const auto &t : o.telescopes)
      run.telescopes.push_back({++id, {t.x / 100.0, t.y / 100.0, t.z / 100.0}, t.r / 100.0});
  }
  return run;
}

// Emissions of case `index` (one source). Original count: floor(N / bunch) bunches plus
// one with the fractional probability. Beyond 1.02 x the fiducial radius an isotropic
// source emits only into the cone that contains the fiducial sphere around the origin
// (LightSource::Emit, `frac`), which is a SphereTarget here.
[[nodiscard]] inline std::vector<SourceEmission>
xyzls_emissions(const XyzlsOptions &o, const SkyRun &run, std::size_t index) {
  if (!(o.direction.x != 0.0 || o.direction.y != 0.0 || o.direction.z != 0.0))
    throw std::invalid_argument("xyzls: direction vector must not be zero");
  const double bunch = std::max(1.0, o.bunch_size);
  const double photons = std::max(0.0, o.photons.at(index));
  const Vec3 pos{o.x_cm / 100.0, o.y_cm / 100.0, o.z_cm / 100.0};
  Source s;
  s.frame = gstrickd::Frame::sky;
  s.geometry = Emitter{pos, normalised(o.direction), o.angular};
  s.spectrum = o.spectrum;
  s.pulse = o.pulse;
  s.photons = photons / bunch;
  const double rmax = fiducial_super_radius(run.telescopes);
  if (norm(pos) > 1.02 * rmax) {
    if (std::holds_alternative<Isotropic>(o.angular))
      s.target = SphereTarget{Vec3{}, 1.02 * rmax};
    else if (!std::holds_alternative<FixedAngle>(o.angular))
      throw std::invalid_argument(
          "xyzls: a non-isotropic profile beyond 1.02 x the fiducial radius is not ported "
          "(the original restricts its directions per profile); use isotropic or Cone/Parallel");
  }
  return {SourceEmission{s, o.seed + index, bunch, 0.0}};
}

[[nodiscard]] inline SkyPlan<SourceEmission> plan_xyzls(const XyzlsOptions &o) {
  SkyPlan<SourceEmission> plan;
  plan.run = xyzls_run(o);
  plan.events_per_case = o.events;
  for (std::size_t c = 0; c < o.photons.size(); ++c)
    plan.cases.push_back(xyzls_emissions(o, plan.run, c));
  return plan;
}

inline void write_xyzls(std::ostream &out, const XyzlsOptions &o) {
  write_plan(out, "dlamb-xyzls", plan_xyzls(o));
}

} // namespace dlamb::programs
