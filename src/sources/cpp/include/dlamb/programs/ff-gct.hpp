#pragma once

// Port of sim_telarray/LightEmission/ff-gct.cc: four flashers at the corners of the
// camera of a dual-mirror telescope (GCT/CHEC), firing down onto the secondary mirror.
//
// Geometry: origin at the secondary-mirror centre (telescope centre, as in the original),
// z along the optical axis towards the sky. The flashers sit at (+-xy, +-xy,
// mirror_camera_distance + depth) and point downwards (-90 deg + inclination, radially
// outwards). Their light must reach the secondary, so the photons are written in the
// telescope frame T at their emission points; the optics stage applies the secondary
// reflection. The original keeps the flasher photons inside the telescope in the same way
// (Bypass_Optics=1: the primary is bypassed, the secondary is not).

#include "dlamb/programs/stream_output.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace dlamb::programs::ff_gct {

// Options with the names of ff-gct.cc. Defaults are the values in the code.
struct Options {
  int events{1};                      // --events, ff-gct.cc:129, 181-182
  std::vector<double> photons{2.5e6}; // --photons [4 pi per flasher], ff-gct.cc:120, 159-160
  double bunch_size{1.0};             // --bunchsize, ff-gct.cc:122, 161-162
  double flasher_xy_cm{14.62};        // --flasher-xy [cm], ff-gct.cc:123, 163-164
  double flasher_depth_cm{2.46};      // --flasher-depth [cm], ff-gct.cc:124, 165-166
  double inclination_deg{0.0};        // --flasher-inclination [deg], ff-gct.cc:125, 167-168
  double mirror_camera_distance_cm{356.1075 - 305.03801}; // ff-gct.cc:126, 169-170
  std::string spectrum{"400"};                            // --spectrum, ff-gct.cc:127, 171-172
  std::string lightpulse{"Simple:0"};                     // --lightpulse, ff-gct.cc:127, 173-174
  std::string angular{"isotropic"}; // --angular-distribution, ff-gct.cc:127, 175-176
  std::array<double, 5> delays{};   // --delays [ns], ff-gct.cc:130, 183-185
  std::string fire{"all"};          // --fire all|alternate|demo|1..4, ff-gct.cc:131, 186-208
  int run_number{1};                // --run, ff-gct.cc:129, 177-178
  double altitude_m{2100.0};        // --altitude [m], ff-gct.cc:118, 157-158
  int atmosphere{1};                // --atmosphere, ff-gct.cc:128, 179-180 (numbers only)
  std::uint64_t seed{1};            // new: Philox seed
};

// Secondary fiducial radius [cm] (ff-gct.cc:119, not an option of the original).
inline constexpr double secondary_fiducial_cm = 120.0;

inline void validate(const Options &o) {
  if (o.events < 0) {
    bad("ff-gct: --events must be >= 0");
  }
  if (o.photons.empty()) {
    bad("ff-gct: --photons needs at least one value");
  }
  for (const double p : o.photons) {
    if (!(p >= 0.0)) {
      bad("ff-gct: --photons must be >= 0");
    }
  }
  if (!(o.flasher_xy_cm >= 0.0) || !(o.mirror_camera_distance_cm > 0.0)) {
    bad("ff-gct: positions must be >= 0 (xy) and > 0 (mirror distance)");
  }
  if (o.fire != "all" && o.fire != "alternate" && o.fire != "demo" &&
      !(o.fire.size() == 1 && o.fire[0] >= '1' && o.fire[0] <= '4')) {
    bad("ff-gct: --fire must be all, alternate, demo or 1..4");
  }
}

// Flasher i (ls1..ls4 of ff-gct.cc:218-229) in the T frame, in metres.
inline Emitter flasher(const Options &o, int i) {
  const double sx = i < 2 ? 1.0 : -1.0;        // ls1, ls2 at +x; ls3, ls4 at -x
  const double sy = (i % 2 == 0) ? 1.0 : -1.0; // ls1, ls3 at +y; ls2, ls4 at -y
  const double xy = o.flasher_xy_cm;
  const double x = sx * xy;
  const double y = sy * xy;
  const double z = o.mirror_camera_distance_cm + o.flasher_depth_cm;
  // DirectionVect(alt, az) of ff-gct.cc:220: alt = -90 deg + inclination, az = atan2(y, x).
  constexpr double deg = 3.141592653589793238462643383279502884 / 180.0;
  const double alt = (-90.0 + o.inclination_deg) * deg;
  const double az = std::atan2(y, x);
  const Vec3 pointing{std::cos(alt) * std::cos(az), std::cos(alt) * std::sin(az), std::sin(alt)};
  return Emitter{Vec3{0.01 * x, 0.01 * y, 0.01 * z}, pointing, parse_angular(o.angular)};
}

// Sources of all events (ff-gct.cc:233-263). Each emission is one flasher firing with
// time offset delays[k] (the `when` of LightSource::Emit).
inline std::vector<Event> build_sources(const Options &o) {
  validate(o);
  const double bunch = o.bunch_size < 1.0 ? 1.0 : o.bunch_size; // IactLightEmission.cc:1514
  const double radius_m = 0.01 * secondary_fiducial_cm;
  const Angular angular = parse_angular(o.angular);
  const Spectrum spectrum = parse_spectrum(o.spectrum);
  const Pulse pulse = parse_pulse(o.lightpulse);

  // "demo" without delays uses the decent demo values (ff-gct.cc:197-206).
  std::array<double, 5> delays = o.delays;
  if (o.fire == "demo" &&
      std::all_of(delays.begin(), delays.end(), [](double d) { return d == 0.0; })) {
    delays = {0.0, 10.0, 20.0, 30.0, 50.0};
  }

  std::vector<Event> result;
  std::uint64_t emission_index = 0;
  for (const auto &photons : o.photons) {
    auto emit_flasher = [&](Event &event, int i, double delay) {
      const Emitter emitter = flasher(o, i);
      Source source;
      source.frame = gstrickd::Frame::telescope;
      source.geometry = emitter;
      source.spectrum = spectrum;
      source.pulse = pulse;
      source.photons = photons / bunch;
      source.target = fiducial_target(emitter.position, radius_m, angular);
      event.emissions.push_back({source, derive_seed(o.seed, emission_index++), bunch, delay});
    };
    for (int iev = 0; iev < o.events; ++iev) {
      Event event;
      event.number = static_cast<std::uint32_t>(result.size() + 1);
      if (o.fire == "all") {
        for (int i = 0; i < 4; ++i) {
          emit_flasher(event, i, delays[static_cast<std::size_t>(i)]);
        }
      } else if (o.fire == "alternate") {
        const int i = iev % 4;
        emit_flasher(event, i, delays[static_cast<std::size_t>(i)]);
      } else if (o.fire == "demo") {
        for (int i = 0; i < 4; ++i) {
          emit_flasher(event, i, delays[static_cast<std::size_t>(i)]);
        }
        for (int i = 0; i < 4; ++i) {
          emit_flasher(event, i, delays[4]);
        }
      } else {
        const int i = o.fire[0] - '1';
        emit_flasher(event, i, delays[static_cast<std::size_t>(i)]);
      }
      result.push_back(std::move(event));
    }
  }
  return result;
}

inline RunSetup run_setup(const Options &o) {
  RunSetup run;
  run.producer = "dlamb-ff-gct";
  run.run_number = o.run_number;
  run.observation_level_m = o.altitude_m;
  run.atmosphere = o.atmosphere;
  run.telescope_radius_m = 0.01 * secondary_fiducial_cm;
  run.frame = gstrickd::Frame::telescope;
  return run;
}

inline void write_stream(std::ostream &out, const Options &o) {
  dlamb::programs::write_stream(out, run_setup(o), build_sources(o));
}

} // namespace dlamb::programs::ff_gct
