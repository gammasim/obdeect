#pragma once

// Port of LightEmission fake-muon.cc: Cherenkov light of a single muon-like track, emitted
// at sky level (frame sky, altitude above the observation level) and filtered to a
// telescope sphere at the observation-level origin. Refraction and extinction are applied
// downstream (dlufd), not here.

#include "dlamb/cli.hpp"
#include "dlamb/emission_stream.hpp"
#include "dlamb/track.hpp"

#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dlamb::fake_muon {

// Refractive index minus one as a function of altitude above the observation level [m]
// (LightEmission: run.RefIdx(dz)). Linear between nodes, constant outside.
struct IndexProfile {
  std::vector<double> altitude_m;
  std::vector<double> index_minus_one;
};

// Options of fake-muon.cc (line numbers in the original). Lengths in metres, angles in
// degrees, wavelengths in nanometres.
struct Options {
  int run_number{1};                  // --run (fake-muon.cc:260)
  std::uint64_t seed{1};              // --seed (new: explicit random seed)
  std::uint32_t events{1};            // -N, --events (fake-muon.cc:257)
  double observation_level_m{2150.0}; // -a (fake-muon.cc:239)
  double telescope_radius_m{15.0};    // -r (fake-muon.cc:240)
  double impact_x_m{0.0};             // -x (fake-muon.cc:241)
  double impact_y_m{0.0};             // -y (fake-muon.cc:242)
  double rscat_m{0.0};                // --rscat (fake-muon.cc:243); must be 0
  int nscat{1};                       // --nscat (fake-muon.cc:244); must be 1
  double h_start_m{3000.0};           // -h (fake-muon.cc:245)
  double theta_deg{0.0};              // --theta (fake-muon.cc:249)
  double phi_deg{0.0};                // --phi (fake-muon.cc:249)
  double index_minus_one{0.000292};   // -n (fake-muon.cc:250)
  IndexProfile index_profile{};       // -A (fake-muon.cc:258); empty = constant
  double bunch_size{5.0};             // -b (fake-muon.cc:247); no effect here
  double wavelength_min_nm{240.0};    // -l (fake-muon.cc:248)
  double wavelength_max_nm{1000.0};   // -L (fake-muon.cc:248)
  double charge{-1.0};                // -q (fake-muon.cc:251)
  double beta{1.0};                   // --beta (fake-muon.cc:252)
  double mass_gev{0.105658};          // -m (fake-muon.cc:253)
  int primary_id{5};                  // --primary (fake-muon.cc:262)
  bool multiple_scattering{false};    // --enable-scatter (fake-muon.cc:259)
};

// Fine-structure constant (fake-muon.cc:256).
inline constexpr double alpha_fine = 0.0072973526;

// Mass [GeV] and charge of a CORSIKA particle type (port of get_mass_charge,
// fake-muon.cc:63-179). Unknown types give mass 1.2 GeV and charge 0.
inline void mass_charge(int type, double &mass, double &charge) {
  charge = mass = 0.0;
  if ((type > 1 && type < 100) || (type > 110 && type < 199)) {
    switch (type) {
    case 3:
      charge = -1;
      mass = 0.000511;
      break; // electron
    case 2:
      charge = 1;
      mass = 0.000511;
      break; // positron
    case 6:
      charge = -1;
      mass = 0.105658;
      break; // mu-
    case 5:
      charge = 1;
      mass = 0.105658;
      break; // mu+
    case 8:
      charge = 1;
      mass = 0.13957;
      break; // pi+
    case 11:
      charge = 1;
      mass = 0.493677;
      break; // K+
    case 14:
      charge = 1;
      mass = 0.938272;
      break; // proton
    case 13:
    case 25:
      charge = 0;
      mass = 0.939565;
      break;
    case 9:
      charge = -1;
      mass = 0.13957;
      break; // pi-
    case 12:
      charge = -1;
      mass = 0.493677;
      break; // K-
    case 10:
    case 16:
      charge = 0;
      mass = 0.497648;
      break;
    case 15:
      charge = -1;
      mass = 0.938272;
      break; // antiproton
    case 66:
    case 67:
    case 68:
    case 69:
    case 133:
    case 134:
      charge = 0;
      mass = 1e-11;
      break; // neutrinos
    case 130:
      charge = 1;
      mass = 1.777;
      break; // tau+
    case 131:
      charge = -1;
      mass = 1.777;
      break; // tau-
    case 7:
      charge = 0;
      mass = 0.134;
      break; // pi0
    case 17:
      charge = 0;
      mass = 0.54751;
      break; // eta
    case 50:
      charge = 0;
      mass = 0.78265;
      break; // omega(782)
    case 51:
      charge = 0;
      mass = 0.7755;
      break;
    case 52:
      charge = 1;
      mass = 0.7755;
      break;
    case 53:
      charge = -1;
      mass = 0.7755;
      break;
    default:
      charge = 0;
      mass = 1.2;
      break; // fall-back
    }
  } else if (type >= 100) {
    charge = type % 100;
    mass = (0.931494 - 0.000511) * static_cast<int>(type / 100);
    if (type == 402)
      mass = 0.931845 * 4.0; // 4He
  }
}

// Lorentz factor, momentum and energy of the track (fake-muon.cc:404-410). Energy [GeV] is
// 1e-9 for a photon-like track (beta 1 or no mass).
[[nodiscard]] inline double energy_gev(const Options &o) {
  if (o.beta < 1.0 && o.mass_gev > 0.0) {
    const double gamma = 1.0 / std::sqrt(1.0 - o.beta * o.beta);
    return gamma * o.mass_gev;
  }
  return 1e-9;
}

// Throws std::invalid_argument for options the port cannot or should not run.
inline void validate(const Options &o) {
  const auto fail = [](const std::string &message) { throw std::invalid_argument(message); };
  if (o.multiple_scattering)
    fail("multiple scattering is not ported (needs the density profile; fake-muon.cc:534-589)");
  if (o.rscat_m != 0.0)
    fail("--rscat is not supported: the impact randomisation is disabled in fake-muon.cc:468");
  if (o.nscat != 1)
    fail("--nscat must be 1 (each event is used once)");
  if (o.events < 1)
    fail("events must be >= 1");
  if (!(o.telescope_radius_m > 0.0) || !std::isfinite(o.telescope_radius_m))
    fail("telescope radius must be positive");
  if (!std::isfinite(o.impact_x_m) || !std::isfinite(o.impact_y_m))
    fail("impact position must be finite");
  if (!(o.h_start_m > 0.0) || !std::isfinite(o.h_start_m))
    fail("starting altitude must be positive");
  if (!(o.theta_deg >= 0.0 && o.theta_deg < 90.0))
    fail("theta must lie in [0, 90) deg");
  if (!std::isfinite(o.phi_deg))
    fail("phi must be finite");
  if (!(o.index_minus_one > 0.0 && o.index_minus_one < 0.0003))
    fail("index of refraction minus one must lie in (0, 0.0003)");
  if (!(o.beta > 0.0 && o.beta <= 1.0))
    fail("beta must lie in (0, 1]");
  if (!(o.mass_gev >= 0.0005 && o.mass_gev <= 100.0))
    fail("mass must lie in [0.0005, 100] GeV/c^2");
  if (o.charge == 0.0)
    fail("a neutral particle emits no Cherenkov light");
  if (!(o.wavelength_min_nm >= 200.0 && o.wavelength_min_nm <= o.wavelength_max_nm &&
        o.wavelength_max_nm <= 2000.0))
    fail("wavelengths must satisfy 200 <= wlmin <= wlmax <= 2000 nm");
  if (!(o.bunch_size > 0.0))
    fail("bunch size must be positive");
  const auto &p = o.index_profile;
  if (!p.altitude_m.empty()) {
    if (p.altitude_m.size() < 2 || p.altitude_m.size() != p.index_minus_one.size())
      fail("an index profile needs at least two nodes and one value per node");
    for (std::size_t i = 0; i < p.altitude_m.size(); ++i) {
      if (!std::isfinite(p.altitude_m[i]) || !std::isfinite(p.index_minus_one[i]) ||
          p.index_minus_one[i] <= 0.0 || (i > 0 && !(p.altitude_m[i] > p.altitude_m[i - 1])))
        fail("index profile altitudes must increase and index values be positive");
    }
  }
}

// Index minus one at altitude h above the observation level.
[[nodiscard]] inline double index_at(const Options &o, double h) {
  const auto &p = o.index_profile;
  if (p.altitude_m.empty())
    return o.index_minus_one;
  if (h <= p.altitude_m.front())
    return p.index_minus_one.front();
  if (h >= p.altitude_m.back())
    return p.index_minus_one.back();
  std::size_t i = 1;
  while (p.altitude_m[i] < h)
    ++i;
  const double f = (h - p.altitude_m[i - 1]) / (p.altitude_m[i] - p.altitude_m[i - 1]);
  return p.index_minus_one[i - 1] + f * (p.index_minus_one[i] - p.index_minus_one[i - 1]);
}

// Photons per metre of track emitted into 4 pi (Frank-Tamm; fake-muon.cc:504 and 594 with
// the cm to m factor of 100 and the nm/cm factor of 1e7).
[[nodiscard]] inline double photons_per_metre(double charge, double sin2_theta, double wl1,
                                              double wl2) {
  return 2.0 * pi * alpha_fine * charge * charge * sin2_theta * (1.0 / wl1 - 1.0 / wl2) * 1e9;
}

// Track segments of one event (the impact position is the same for all events). Without an
// index profile the track is one segment with a constant index (LightEmission atmosphere
// "const", fake-muon.cc:486). With a profile it is cut into 10 m steps, each with the index
// at its mid-point (fake-muon.cc:529-611). A step with cos(theta) > 1 emits nothing; the
// original skips it without advancing the track (a bug, not reproduced: the geometry stays
// consistent here).
[[nodiscard]] inline std::vector<Track> tracks(const Options &o) {
  validate(o);
  const double theta = o.theta_deg * pi / 180.0, phi = o.phi_deg * pi / 180.0;
  const double sin_t = std::sin(theta), cos_t = std::cos(theta);
  const Vec3 direction{-sin_t * std::cos(phi), -sin_t * std::sin(phi), -cos_t};
  const double dist = o.h_start_m / cos_t;
  const double c = speed_of_light_m_per_ns;
  const Spectrum spectrum = o.wavelength_min_nm == o.wavelength_max_nm
                                ? Spectrum{Fixed{o.wavelength_min_nm}}
                                : Spectrum{InverseSquare{o.wavelength_min_nm, o.wavelength_max_nm}};
  const auto make = [&](Vec3 start, double length, double cos_ct, double time) {
    const double sin2 = 1.0 - cos_ct * cos_ct;
    Track t;
    t.start = start;
    t.direction = direction;
    t.length_m = length;
    t.cos_theta = cos_ct;
    t.beta = o.beta;
    t.time_ns = time;
    t.spectrum = spectrum;
    t.photons =
        photons_per_metre(o.charge, sin2, o.wavelength_min_nm, o.wavelength_max_nm) * length;
    t.target = SphereTarget{{0.0, 0.0, 0.0}, o.telescope_radius_m};
    return t;
  };
  std::vector<Track> result;
  const Vec3 impact{o.impact_x_m, o.impact_y_m, 0.0};
  if (o.index_profile.altitude_m.empty()) {
    const double n = 1.0 + o.index_minus_one;
    const double cos_ct = 1.0 / (n * o.beta);
    if (cos_ct > 1.0)
      return result;
    const Vec3 start = impact - direction * dist; // upper end, at altitude h_start
    result.push_back(make(start, dist, cos_ct, -dist / (o.beta * c)));
    return result;
  }
  const auto steps = static_cast<std::size_t>(std::ceil(dist / 10.0));
  const double step = dist / static_cast<double>(steps);
  Vec3 lower = impact;
  double travelled = 0.0;
  for (std::size_t k = 0; k < steps; ++k) {
    const Vec3 upper = lower - direction * step;
    const Vec3 mid = lower - direction * (0.5 * step);
    travelled += step;
    const double n = 1.0 + index_at(o, mid.z);
    const double cos_ct = 1.0 / (n * o.beta);
    if (cos_ct <= 1.0)
      result.push_back(make(upper, step, cos_ct, -travelled / (o.beta * c)));
    lower = upper;
  }
  return result;
}

// Writes the complete stream of fake-muon. Source k of the run uses seed source_seed(seed, k).
// Returns the expected and written photon weights (diagnostics).
inline StreamSummary write(std::ostream &out, const Options &o,
                           const std::string &producer = "dlamb-fake-muon") {
  const auto segments = tracks(o);
  const TelescopeSphere telescope{1, {0.0, 0.0, 0.0}, o.telescope_radius_m};
  RunInfo run;
  run.run_number = o.run_number;
  run.observation_level_m = o.observation_level_m;
  run.wavelength_min_nm = o.wavelength_min_nm;
  run.wavelength_max_nm = o.wavelength_max_nm;
  run.telescopes = {telescope};
  RunWriter writer(out, producer, run);
  double expected_per_event = 0.0;
  for (const auto &t : segments)
    expected_per_event += expected_track_weight(t);
  for (std::uint32_t event = 1; event <= o.events; ++event) {
    writer.begin_event({event, o.primary_id, energy_gev(o), o.impact_x_m, o.impact_y_m});
    writer.begin_telescope(telescope, gstrickd::Frame::sky);
    std::vector<Photon> photons;
    for (std::size_t k = 0; k < segments.size(); ++k) {
      photons.clear();
      emit_track(segments[k], source_seed(o.seed, k), event, telescope.id, photons);
      for (const auto &p : photons)
        writer.push(p);
    }
    writer.end_telescope();
    writer.end_event();
  }
  writer.finish();
  return {writer.photons(), writer.weight(), expected_per_event * o.events};
}

} // namespace dlamb::fake_muon
