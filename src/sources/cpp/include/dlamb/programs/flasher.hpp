#pragma once

// Flasher photon streams of two types, selected by Type. Replace the earlier ff-1m (dish)
// and ff-gct (secondary) programs; options keep their names.
//
// dish: flat-fielding flasher on the axis of a single-reflector telescope, in
// front of the camera. The flasher sits at (x, y, distance) above the telescope centre, the
// origin at observation level. Photons are written in the sky frame at their emission
// points (or, with frame camera, on the camera plane); the atmosphere stage moves the sky
// frame photons to the ground.
//
// secondary: four flashers at the corners of the camera of a dual-mirror
// telescope (GCT/CHEC), firing down onto the secondary mirror. Origin at the secondary
// centre, z along the optical axis towards the sky. The flashers sit at
// (+-xy, +-xy, mirror_camera_distance + depth) and point downwards (-90 deg + inclination,
// radially outwards). Their light must reach the secondary, so the photons are written in
// the telescope frame T at their emission points.

#include "dlamb/programs/stream_output.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

namespace dlamb::programs::flasher {

enum class Type { dish, secondary };

// --type dish|secondary.
inline Type parse_type(const std::string &text) {
  if (text == "dish") {
    return Type::dish;
  }
  if (text == "secondary") {
    return Type::secondary;
  }
  bad("flasher: --type must be dish or secondary, not '" + text + "'");
}

// Rejects an option of `owner` used with the other type.
inline void require_type(Type type, Type owner, const std::string &option) {
  if (type != owner) {
    bad(option + " applies only to --type " + (owner == Type::dish ? "dish" : "secondary"));
  }
}

// Option names as in the earlier ff-1m (dish) and ff-gct (secondary) programs.
struct Options {
  Type type{Type::dish};
  // dish only: sky (emission points) or camera (photons on the camera plane, no atmosphere).
  gstrickd::Frame frame{gstrickd::Frame::sky};
  // Both types.
  std::vector<double> photons{2.5e6}; // --photons [4 pi per flasher]
  double bunch_size{1.0};             // --bunchsize [photons]
  std::string spectrum{"400"};        // --spectrum [nm]
  std::string lightpulse{"Simple:0"}; // --lightpulse
  std::string angular{"isotropic"};   // --angular-distribution
  int run_number{1};                  // --run
  std::optional<double> altitude_m;   // --altitude [m]: 2150 (dish), 2100 (secondary)
  int atmosphere{1};                  // --atmosphere (numbers only)
  std::uint64_t seed{1};              // new: Philox seed (counter-based random numbers)
  // dish: per photon level (a short list repeats its last value).
  // secondary: one value.
  std::vector<int> events{1}; // --events
  // dish only.
  double x_cm{0.0};               // --xy x,y (first value) [cm]
  double y_cm{0.0};               // --xy x,y (second value) [cm]
  double distance_cm{1600.0};     // --distance [cm]
  double camera_radius_cm{120.0}; // --camera-radius [cm]
  // secondary only.
  double flasher_xy_cm{14.62};                            // --flasher-xy [cm]
  double flasher_depth_cm{2.46};                          // --flasher-depth [cm]
  double inclination_deg{0.0};                            // --flasher-inclination [deg]
  double mirror_camera_distance_cm{356.1075 - 305.03801}; // [cm]
  std::array<double, 5> delays{};                         // --delays [ns]
  std::string fire{"all"};                                // --fire all|alternate|demo|1..4
};

// Secondary fiducial radius [cm] (not an option).
inline constexpr double secondary_fiducial_cm = 120.0;

inline double observation_level_m(const Options &o) {
  if (o.altitude_m) {
    return *o.altitude_m;
  }
  return o.type == Type::dish ? 2150.0 : 2100.0;
}

inline void validate(const Options &o) {
  if (o.events.empty()) {
    bad("flasher: --events needs at least one value");
  }
  for (const int n : o.events) {
    if (n < 0) {
      bad("flasher: --events must be >= 0");
    }
  }
  if (o.photons.empty()) {
    bad("flasher: --photons needs at least one value");
  }
  for (const double p : o.photons) {
    if (!(p >= 0.0) || !std::isfinite(p)) {
      bad("flasher: --photons must be finite and >= 0");
    }
  }
  // Bunch sizes below 1 are raised to 1 when the sources are built.
  if (!std::isfinite(o.bunch_size)) {
    bad("flasher: --bunchsize must be finite");
  }
  if (o.frame != gstrickd::Frame::sky && o.frame != gstrickd::Frame::camera) {
    bad("flasher: --frame must be sky or camera");
  }
  if (o.type == Type::secondary && o.frame != gstrickd::Frame::sky) {
    bad("flasher: --frame applies only to --type dish");
  }
  if (o.type == Type::dish) {
    if (!(o.distance_cm > 0.0) || !(o.camera_radius_cm > 0.0)) {
      bad("flasher: --distance and --camera-radius must be > 0 (--type dish)");
    }
  } else {
    if (o.events.size() != 1) {
      bad("flasher: --events takes one value with --type secondary");
    }
    if (!(o.flasher_xy_cm >= 0.0) || !(o.mirror_camera_distance_cm > 0.0)) {
      bad("flasher: positions must be >= 0 (xy) and > 0 (mirror distance) (--type secondary)");
    }
    if (o.fire != "all" && o.fire != "alternate" && o.fire != "demo" &&
        !(o.fire.size() == 1 && o.fire[0] >= '1' && o.fire[0] <= '4')) {
      bad("flasher: --fire must be all, alternate, demo or 1..4");
    }
  }
}

// Dish placement. Event i fires one flasher emission per photon level, for
// events[k] events at photon level k. With frame camera the photons are propagated in a
// straight line to the camera plane z = 0 and the time includes path / c.
inline std::vector<Event> dish_sources(const Options &o) {
  // At least one photon per bunch.
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
      source.frame = o.frame;
      source.geometry = Emitter{position, Vec3{0.0, 0.0, -1.0}, angular};
      source.spectrum = parse_spectrum(o.spectrum);
      source.pulse = parse_pulse(o.lightpulse);
      source.photons = o.photons[level] / bunch;
      source.target = target;
      Event event;
      event.number = static_cast<std::uint32_t>(result.size() + 1);
      event.emissions.push_back({source, derive_seed(o.seed, emission++), bunch, 0.0,
                                 o.frame == gstrickd::Frame::camera});
      result.push_back(std::move(event));
    }
  }
  return result;
}

// Flasher i (1..4) in the T frame, in metres.
inline Emitter secondary_emitter(const Options &o, int i) {
  const double sx = i < 2 ? 1.0 : -1.0;        // ls1, ls2 at +x; ls3, ls4 at -x
  const double sy = (i % 2 == 0) ? 1.0 : -1.0; // ls1, ls3 at +y; ls2, ls4 at -y
  const double xy = o.flasher_xy_cm;
  const double x = sx * xy;
  const double y = sy * xy;
  const double z = o.mirror_camera_distance_cm + o.flasher_depth_cm;
  // Pointing: alt = -90 deg + inclination, az = atan2(y, x).
  constexpr double deg = 3.141592653589793238462643383279502884 / 180.0;
  const double alt = (-90.0 + o.inclination_deg) * deg;
  const double az = std::atan2(y, x);
  const Vec3 pointing{std::cos(alt) * std::cos(az), std::cos(alt) * std::sin(az), std::sin(alt)};
  return Emitter{Vec3{0.01 * x, 0.01 * y, 0.01 * z}, pointing, parse_angular(o.angular)};
}

// Secondary placement. Each emission is one flasher firing with time
// offset delays[k] (time offset of the emission).
inline std::vector<Event> secondary_sources(const Options &o) {
  const double bunch = o.bunch_size < 1.0 ? 1.0 : o.bunch_size;
  const double radius_m = 0.01 * secondary_fiducial_cm;
  const Angular angular = parse_angular(o.angular);
  const Spectrum spectrum = parse_spectrum(o.spectrum);
  const Pulse pulse = parse_pulse(o.lightpulse);
  const int events = o.events.front();

  // "demo" without delays uses the decent demo values.
  std::array<double, 5> delays = o.delays;
  if (o.fire == "demo" &&
      std::all_of(delays.begin(), delays.end(), [](double d) { return d == 0.0; })) {
    delays = {0.0, 10.0, 20.0, 30.0, 50.0};
  }

  std::vector<Event> result;
  std::uint64_t emission_index = 0;
  for (const auto &photons : o.photons) {
    auto emit_flasher = [&](Event &event, int i, double delay) {
      const Emitter emitter = secondary_emitter(o, i);
      Source source;
      source.frame = gstrickd::Frame::telescope;
      source.geometry = emitter;
      source.spectrum = spectrum;
      source.pulse = pulse;
      source.photons = photons / bunch;
      source.target = fiducial_target(emitter.position, radius_m, angular);
      event.emissions.push_back({source, derive_seed(o.seed, emission_index++), bunch, delay});
    };
    for (int iev = 0; iev < events; ++iev) {
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

// Sources of all events of the selected type.
inline std::vector<Event> build_sources(const Options &o) {
  validate(o);
  return o.type == Type::dish ? dish_sources(o) : secondary_sources(o);
}

inline RunSetup run_setup(const Options &o) {
  RunSetup run;
  run.producer = "dlamb-flasher";
  run.run_number = o.run_number;
  run.observation_level_m = observation_level_m(o);
  run.atmosphere = o.atmosphere;
  if (o.type == Type::dish) {
    run.telescope_radius_m = 0.01 * o.camera_radius_cm;
    run.frame = o.frame;
  } else {
    run.telescope_radius_m = 0.01 * secondary_fiducial_cm;
    run.frame = gstrickd::Frame::telescope;
  }
  return run;
}

// Writes the complete gstrickd stream for the options.
inline void write_stream(std::ostream &out, const Options &o) {
  dlamb::programs::write_stream(out, run_setup(o), build_sources(o));
}

} // namespace dlamb::programs::flasher
