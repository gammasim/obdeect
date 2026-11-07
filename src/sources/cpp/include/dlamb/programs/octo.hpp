#pragma once

// Port of LightEmission octo.cc: an isotropic point light moved along a figure-eight at
// 1000 m above the observation level, emitting at each position. Photons are written at
// sky level; only those that reach the telescope sphere (at the observation-level origin)
// are generated (sphere target, LightEmission keeps the same photons as bunches).

#include "dlamb/cli.hpp"
#include "dlamb/emission_stream.hpp"
#include "dlamb/source.hpp"

#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dlamb::octo {

// Options of octo.cc. octo.cc has no command-line options; the CLI maps the fields below.
struct Options {
  double observation_level_m{2150.0}; // alt_m (octo.cc:76)
  double telescope_radius_m{12.5};    // r_tel_m (octo.cc:77)
  double distance_m{1000.0};          // dist (octo.cc:78): light altitude above obs. level
  double photons_per_position{8e7};   // nphot_per_pos (octo.cc:79), photons into 4 pi
  double ring_radius_m{20.0};         // rc (octo.cc:86): radius of the two circles
  double wavelength_nm{400.0};        // SimpleLightSource default (IactLightEmission.hh:591)
  std::uint64_t seed{1};              // new: explicit random seed
  int run_number{1};                  // Run(1, ...) (octo.cc:83)
};

// One emission position: light source at `position` emitting at `time_ns`.
struct Position {
  Vec3 position{};
  double time_ns{};
};

// Positions of octo.cc:107-124 in the same order, with the same loop variables (41 positions
// on the first circle at t = 0 .. 20 ns, 40 on the second at t = 0.5 .. 20 ns, shifted by
// 20 ns). The circle makes one turn in 20 ns.
[[nodiscard]] inline std::vector<Position> figure_eight(const Options &o) {
  const double rc = o.ring_radius_m, d = o.distance_m;
  std::vector<Position> result;
  for (double t = 0.0; t < 20.1; t += 0.5)
    result.push_back(
        {{rc - rc * std::cos(2.0 * pi * t / 20.0), rc * std::sin(2.0 * pi * t / 20.0), d}, t});
  for (double t = 0.5; t < 20.1; t += 0.5)
    result.push_back(
        {{-rc + rc * std::cos(2.0 * pi * t / 20.0), rc * std::sin(2.0 * pi * t / 20.0), d},
         t + 20.0});
  return result;
}

// Throws std::invalid_argument for invalid options.
inline void validate(const Options &o) {
  const auto fail = [](const std::string &message) { throw std::invalid_argument(message); };
  if (!(o.photons_per_position >= 0.0) || !std::isfinite(o.photons_per_position))
    fail("photons per position must be finite and >= 0");
  if (!(o.telescope_radius_m > 0.0) || !std::isfinite(o.telescope_radius_m))
    fail("telescope radius must be positive");
  if (!(o.distance_m > o.telescope_radius_m) || !std::isfinite(o.distance_m))
    fail("the light must be above the telescope sphere");
  if (!(o.ring_radius_m >= 0.0) || !std::isfinite(o.ring_radius_m))
    fail("ring radius must be finite and >= 0");
  if (!(o.wavelength_nm > 0.0) || !std::isfinite(o.wavelength_nm))
    fail("wavelength must be positive");
}

// Source of one emission position (isotropic, instantaneous, sphere target).
[[nodiscard]] inline Source source_at(const Options &o, const Position &p) {
  Source s;
  s.frame = gstrickd::Frame::sky;
  s.geometry = Emitter{p.position, {0.0, 0.0, -1.0}, Isotropic{}};
  s.spectrum = Fixed{o.wavelength_nm};
  s.pulse = Fixed{p.time_ns};
  s.photons = o.photons_per_position;
  s.target = SphereTarget{{0.0, 0.0, 0.0}, o.telescope_radius_m};
  return s;
}

// Writes the complete stream of one event with all positions. Position k is source k with
// seed source_seed(seed, k). Returns expected and written photon weights.
inline StreamSummary write(std::ostream &out, const Options &o,
                           const std::string &producer = "dlamb-octo") {
  validate(o);
  const TelescopeSphere telescope{1, {0.0, 0.0, 0.0}, o.telescope_radius_m};
  RunInfo run;
  run.run_number = o.run_number;
  run.observation_level_m = o.observation_level_m;
  run.wavelength_min_nm = o.wavelength_nm;
  run.wavelength_max_nm = o.wavelength_nm;
  run.telescopes = {telescope};
  RunWriter writer(out, producer, run);
  writer.begin_event({1, 1, 1e-9, 0.0, 0.0});
  writer.begin_telescope(telescope, gstrickd::Frame::sky);
  double expected = 0.0;
  const auto positions = figure_eight(o);
  for (std::size_t k = 0; k < positions.size(); ++k) {
    const Source s = source_at(o, positions[k]);
    validate(s);
    const std::uint64_t seed = source_seed(o.seed, k);
    expected += expected_count(s);
    const std::uint64_t n = photon_count(s, seed, 1, telescope.id);
    for (std::uint64_t i = 0; i < n; ++i)
      writer.push(sample_photon(s, seed, 1, telescope.id, i));
  }
  writer.end_telescope();
  writer.end_event();
  writer.finish();
  return {writer.photons(), writer.weight(), expected};
}

} // namespace dlamb::octo
