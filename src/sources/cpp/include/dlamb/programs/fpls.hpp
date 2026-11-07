#pragma once

// Port of sim_telarray/LightEmission/fpls.cc: a focal-plane light source (LED) at
// 0.7 times the camera radius above the telescope centre, tilted 15 deg from the vertical
// and rotated about the vertical axis once per rotation period.
//
// Geometry: origin at the telescope centre on the observation level, z up (sky frame).
// Photons are written at their emission points and directions; the original propagates
// them to the observation plane, which sim_telarray uses as the camera plane for
// Bypass_Optics=2. The atmosphere stage (G, refraction) and a projecting optics stage
// reproduce it.
//
// The original reads its spectrum (ls_spec.lis) and beam profile (ls_beamed.lis) from
// files that are not in this checkout. Here they are options; the defaults are a fixed
// 400 nm and an isotropic beam, so the default output is not the original's.

#include "dlamb/programs/stream_output.hpp"

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace dlamb::programs::fpls {

// Options. Values in comments are the constants of fpls.cc (no options in the original).
struct Options {
  double altitude_m{2100.0};          // fpls.cc:67
  double camera_radius_m{1.4};        // fpls.cc:68 (fiducial sphere, LED height 0.7 R)
  double photons{5000.0};             // fpls.cc:70 [4 pi per time step]
  double bunch_size{1.0};             // fpls.cc:71
  double rotation_period_ns{60.0};    // fpls.cc:72
  double time_step_ns{0.25};          // fpls.cc:87
  std::string spectrum{"400"};        // fpls.cc:82 uses ls_spec.lis (file not available)
  std::string lightpulse{"Simple:0"}; // fpls.cc:82
  std::string angular{"isotropic"};   // fpls.cc:82 uses ls_beamed.lis (file not available)
  int run_number{1};                  // fpls.cc:76
  int atmosphere{1};                  // fpls.cc:73
  std::uint64_t seed{1};              // new: Philox seed
};

// Tilt of the LED axis from the vertical [deg] (fpls.cc:83).
inline constexpr double tilt_deg = 15.0;

inline void validate(const Options &o) {
  if (!(o.camera_radius_m > 0.0) || !(o.photons >= 0.0) || !(o.bunch_size > 0.0) ||
      !(o.rotation_period_ns > 0.0) || !(o.time_step_ns > 0.0)) {
    bad("fpls: radius, photons, bunch size, rotation period and time step must be > 0 "
        "(photons >= 0)");
  }
}

// Sources of the single event: one emission per time step t = 0, step, ... <= period,
// with the LED pointing rotated by 2 pi t / period about the vertical (fpls.cc:86-93).
inline std::vector<Event> build_sources(const Options &o) {
  validate(o);
  constexpr double deg = 3.141592653589793238462643383279502884 / 180.0;
  const double bunch = o.bunch_size < 1.0 ? 1.0 : o.bunch_size; // IactLightEmission.cc:1514
  const Vec3 position{0.0, 0.0, 0.7 * o.camera_radius_m};       // fpls.cc:69
  // DirectionVect(-75 deg, 0): downwards, 15 deg from the vertical (fpls.cc:83).
  const Vec3 axis{std::cos((-90.0 + tilt_deg) * deg), 0.0, std::sin((-90.0 + tilt_deg) * deg)};
  const Angular angular = parse_angular(o.angular);
  const Spectrum spectrum = parse_spectrum(o.spectrum);
  const Pulse pulse = parse_pulse(o.lightpulse);
  const Target target = fiducial_target(position, o.camera_radius_m, angular);

  Event event;
  event.number = 1;
  std::uint64_t emission = 0;
  // Same loop as fpls.cc:87 (the added 0.001 keeps the last step at the period).
  for (double t = 0.0; t < o.rotation_period_ns + 0.001; t += o.time_step_ns) {
    const double phi = t * (2.0 * 3.141592653589793238462643383279502884) / o.rotation_period_ns;
    // Rotation about z (DirectionVect::Rotate_z, fpls.cc:90).
    const Vec3 direction{axis.x * std::cos(phi) - axis.y * std::sin(phi),
                         axis.x * std::sin(phi) + axis.y * std::cos(phi), axis.z};
    Source source;
    source.frame = gstrickd::Frame::sky;
    source.geometry = Emitter{position, direction, angular};
    source.spectrum = spectrum;
    source.pulse = pulse;
    source.photons = o.photons / bunch;
    source.target = target;
    event.emissions.push_back({source, derive_seed(o.seed, emission++), bunch, t});
  }
  return {event};
}

inline RunSetup run_setup(const Options &o) {
  RunSetup run;
  run.producer = "dlamb-fpls";
  run.run_number = o.run_number;
  run.observation_level_m = o.altitude_m;
  run.atmosphere = o.atmosphere;
  run.telescope_radius_m = o.camera_radius_m;
  run.frame = gstrickd::Frame::sky;
  return run;
}

inline void write_stream(std::ostream &out, const Options &o) {
  dlamb::programs::write_stream(out, run_setup(o), build_sources(o));
}

} // namespace dlamb::programs::fpls
