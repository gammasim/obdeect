#pragma once

// Port of sim_telarray/LightEmission/ff-1m.cc: flat-fielding flasher on the axis of a
// single-reflector telescope, in front of the camera.
//
// Geometry: the flasher sits at (x, y, distance) above the telescope centre, which is
// the origin at observation level (Run with one telescope at the origin). Photons are
// written in the sky frame at their emission points. The original propagates them to the
// observation plane (G, with refraction) and sim_telarray then uses that plane as the
// camera plane (Bypass_Optics=2). So the atmosphere stage (G, refraction) and an optics
// stage that only projects (no mirror) reproduce it.

#include "dlamb/programs/stream_output.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace dlamb::programs::camera_flasher {

// Options with the names of ff-1m.cc. Defaults are the values in the code (the
// --help text of the original differs for --altitude: 2100; the code uses 2150).
struct Options {
  std::vector<int> events{1};         // --events, ff-1m.cc:132, 196-200
  std::vector<double> photons{2.5e6}; // --photons [4 pi per flasher], ff-1m.cc:130, 166-167
  double bunch_size{1.0};             // --bunchsize [photons], ff-1m.cc:134, 168-169
  double x_cm{0.0};                   // --xy x,y (first value) [cm], ff-1m.cc:135, 170-179
  double y_cm{0.0};                   // --xy x,y (second value) [cm]
  double distance_cm{1600.0};         // --distance [cm], ff-1m.cc:136, 180-181
  double camera_radius_cm{120.0};     // --camera-radius [cm], ff-1m.cc:129, 182-183
  std::string spectrum{"400"};        // --spectrum [nm], ff-1m.cc:137, 184-185
  std::string lightpulse{"Simple:0"}; // --lightpulse, ff-1m.cc:137, 186-187
  std::string angular{"isotropic"};   // --angular-distribution, ff-1m.cc:137, 188-189
  int run_number{1};                  // --run, ff-1m.cc:138, 190-191
  double altitude_m{2150.0};          // --altitude [m], ff-1m.cc:128, 164-165
  int atmosphere{1};                  // --atmosphere, ff-1m.cc:139, 192-195 (numbers only)
  std::uint64_t seed{1};              // new: Philox seed (the original uses the global RNG)
};

inline void validate(const Options &o) {
  for (const int n : o.events) {
    if (n < 0) {
      bad("camera-flasher: --events must be >= 0");
    }
  }
  if (o.events.empty()) {
    bad("camera-flasher: --events needs at least one value");
  }
  if (o.photons.empty()) {
    bad("camera-flasher: --photons needs at least one value");
  }
  for (const double p : o.photons) {
    if (!(p >= 0.0)) {
      bad("camera-flasher: --photons must be >= 0");
    }
  }
  if (!(o.distance_cm > 0.0) || !(o.camera_radius_cm > 0.0)) {
    bad("camera-flasher: --distance and --camera-radius must be > 0");
  }
}

// Sources of all events. Event i (1-based) fires one flasher emission per photon level,
// for events[k] events at photon level k (the last value repeats, as in ff-1m.cc:218-219).
inline std::vector<Event> build_sources(const Options &o) {
  validate(o);
  // LightSource::Emit uses at least one photon per bunch (IactLightEmission.cc:1514-1515).
  const double bunch = o.bunch_size < 1.0 ? 1.0 : o.bunch_size;
  const Vec3 position{0.01 * o.x_cm, 0.01 * o.y_cm, 0.01 * o.distance_cm};
  const double radius_m = 0.01 * o.camera_radius_cm;
  const Angular angular = parse_angular(o.angular);
  const Target target = fiducial_target(position, radius_m, angular);

  std::vector<Event> result;
  std::uint64_t emission = 0;
  for (std::size_t level = 0; level < o.photons.size(); ++level) {
    const int events = o.events[std::min(level, o.events.size() - 1)];
    for (int i = 0; i < events; ++i) {
      Source source;
      source.frame = gstrickd::Frame::sky;
      source.geometry = Emitter{position, Vec3{0.0, 0.0, -1.0}, angular};
      source.spectrum = parse_spectrum(o.spectrum);
      source.pulse = parse_pulse(o.lightpulse);
      source.photons = o.photons[level] / bunch;
      source.target = target;
      Event event;
      event.number = static_cast<std::uint32_t>(result.size() + 1);
      event.emissions.push_back({source, derive_seed(o.seed, emission++), bunch, 0.0});
      result.push_back(std::move(event));
    }
  }
  return result;
}

inline RunSetup run_setup(const Options &o) {
  RunSetup run;
  run.producer = "dlamb-camera-flasher";
  run.run_number = o.run_number;
  run.observation_level_m = o.altitude_m;
  run.atmosphere = o.atmosphere;
  run.telescope_radius_m = 0.01 * o.camera_radius_cm;
  run.frame = gstrickd::Frame::sky;
  return run;
}

// Writes the complete gstrickd stream for the options.
inline void write_stream(std::ostream &out, const Options &o) {
  dlamb::programs::write_stream(out, run_setup(o), build_sources(o));
}

} // namespace dlamb::programs::camera_flasher
