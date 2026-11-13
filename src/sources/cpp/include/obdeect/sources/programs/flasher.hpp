#pragma once

// Flasher programs: a dish flasher in front of a single-reflector camera, or four flashers of a
// dual-mirror telescope (docs/spec/tool-flasher.md).

#include "obdeect/sources/programs/stream_output.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::sources::programs::flasher {

enum class Type { dish, secondary };

inline Type parse_type(const std::string &text) {
  if (text == "dish")
    return Type::dish;
  if (text == "secondary")
    return Type::secondary;
  throw std::invalid_argument("--type must be dish or secondary, not '" + text + "'");
}

// Rejects an option of one type used with the other.
inline void require_type(Type actual, Type required, const std::string &option) {
  if (actual != required)
    throw std::invalid_argument(option + " applies only to --type " +
                                (required == Type::dish ? "dish" : "secondary"));
}

struct Options {
  Type type{Type::dish};
  obdeect::stream::Frame frame{obdeect::stream::Frame::sky};
  std::vector<double> photons{2.5e6}; // into 4 pi per flasher, one entry per level
  double bunch_size{1.0};
  std::string spectrum{"400"};
  std::string lightpulse{"Simple:0"};
  std::string angular{"isotropic"};
  std::int64_t run_number{1};
  std::optional<double> altitude_m;
  std::int64_t atmosphere{1};
  std::uint64_t seed{1};
  std::vector<int> events{1};
  // dish
  double x_cm{0.0};
  double y_cm{0.0};
  double distance_cm{1600.0};
  double camera_radius_cm{120.0};
  // secondary
  double flasher_xy_cm{14.62};
  double flasher_depth_cm{2.46};
  double inclination_deg{0.0};
  double mirror_camera_distance_cm{356.1075 - 305.03801};
  std::array<double, 5> delays{0.0, 0.0, 0.0, 0.0, 0.0};
  std::string fire{"all"};
};

// Radius of the secondary telescope's fiducial sphere [m].
inline constexpr double secondary_radius_m = 1.2;

inline double observation_level_m(const Options &o) {
  return o.altitude_m ? *o.altitude_m : (o.type == Type::dish ? 2150.0 : 2100.0);
}

namespace detail {

inline Emission make_emission(const Options &o, const Source &source, std::uint64_t index,
                              double bunch, double offset_ns, bool camera_plane) {
  Emission e;
  e.source = source;
  e.seed = derive_seed(o.seed, index);
  e.telescope = 1;
  e.bunch_size = bunch;
  e.time_offset_ns = offset_ns;
  e.to_camera_plane = camera_plane;
  return e;
}

inline std::vector<Event> build_dish(const Options &o, double bunch) {
  if (!(o.distance_cm > 0.0) || !std::isfinite(o.distance_cm))
    throw std::invalid_argument("--distance must be finite and > 0");
  if (!(o.camera_radius_cm > 0.0) || !std::isfinite(o.camera_radius_cm))
    throw std::invalid_argument("--camera-radius must be finite and > 0");
  if (!std::isfinite(o.x_cm) || !std::isfinite(o.y_cm))
    throw std::invalid_argument("--xy must be finite");
  if (o.frame != obdeect::stream::Frame::sky && o.frame != obdeect::stream::Frame::camera)
    throw std::invalid_argument("--frame must be sky or camera");
  if (o.events.empty())
    throw std::invalid_argument("--events needs at least one value");
  for (const int e : o.events)
    if (e < 0)
      throw std::invalid_argument("--events must be >= 0");
  const double radius = o.camera_radius_cm / 100.0;
  const Vec3 position{o.x_cm / 100.0, o.y_cm / 100.0, o.distance_cm / 100.0};
  const Angular angular = parse_angular(o.angular);
  Source source;
  source.frame = o.frame;
  source.geometry = Emitter{position, {0.0, 0.0, -1.0}, angular};
  source.spectrum = parse_spectrum(o.spectrum);
  source.pulse = parse_pulse(o.lightpulse);
  source.target = fiducial_target(position, radius, angular);

  std::vector<Event> events;
  std::uint64_t emission_index = 0;
  for (std::size_t level = 0; level < o.photons.size(); ++level) {
    source.photons = o.photons[level] / bunch;
    validate(source);
    const int count = o.events[std::min(level, o.events.size() - 1)];
    for (int k = 0; k < count; ++k) {
      Event event;
      event.event_number = static_cast<std::int64_t>(events.size()) + 1;
      event.emissions.push_back(make_emission(o, source, emission_index++, bunch, 0.0,
                                              o.frame == obdeect::stream::Frame::camera));
      events.push_back(std::move(event));
    }
  }
  return events;
}

inline std::vector<Event> build_secondary(const Options &o, double bunch) {
  if (o.frame != obdeect::stream::Frame::sky)
    throw std::invalid_argument("--frame applies only to --type dish");
  if (o.events.size() != 1)
    throw std::invalid_argument("--events takes one value with --type secondary");
  if (o.events[0] < 0)
    throw std::invalid_argument("--events must be >= 0");
  if (!(o.flasher_xy_cm >= 0.0) || !std::isfinite(o.flasher_xy_cm) ||
      !std::isfinite(o.flasher_depth_cm) || !std::isfinite(o.inclination_deg))
    throw std::invalid_argument(
        "flasher offset, depth and inclination must be finite, offset >= 0");
  if (!(o.mirror_camera_distance_cm > 0.0) || !std::isfinite(o.mirror_camera_distance_cm))
    throw std::invalid_argument("--mirror-camera-distance must be finite and > 0");
  std::array<double, 5> delays = o.delays;
  for (const double d : delays)
    if (!std::isfinite(d))
      throw std::invalid_argument("--delays must be finite");
  const bool demo = o.fire == "demo";
  if (demo && std::all_of(delays.begin(), delays.end(), [](double d) { return d == 0.0; }))
    delays = {0.0, 10.0, 20.0, 30.0, 50.0};
  int single = -1;
  if (o.fire.size() == 1 && o.fire[0] >= '1' && o.fire[0] <= '4')
    single = o.fire[0] - '1';
  else if (o.fire != "all" && o.fire != "alternate" && !demo)
    throw std::invalid_argument("--fire must be all, alternate, demo or 1..4");

  std::array<Source, 4> sources;
  const double alt = (-90.0 + o.inclination_deg) * pi / 180.0;
  const Spectrum spectrum = parse_spectrum(o.spectrum);
  const Pulse pulse = parse_pulse(o.lightpulse);
  const Angular angular = parse_angular(o.angular);
  for (int i = 0; i < 4; ++i) {
    const double sx = i < 2 ? 1.0 : -1.0, sy = i % 2 == 0 ? 1.0 : -1.0;
    const double x = sx * o.flasher_xy_cm, y = sy * o.flasher_xy_cm;
    const Vec3 position{x / 100.0, y / 100.0,
                        (o.mirror_camera_distance_cm + o.flasher_depth_cm) / 100.0};
    const double az = std::atan2(y, x);
    const Vec3 pointing{std::cos(alt) * std::cos(az), std::cos(alt) * std::sin(az), std::sin(alt)};
    sources[i].frame = obdeect::stream::Frame::telescope;
    sources[i].geometry = Emitter{position, pointing, angular};
    sources[i].spectrum = spectrum;
    sources[i].pulse = pulse;
    sources[i].target = fiducial_target(position, secondary_radius_m, angular);
  }

  std::vector<Event> events;
  std::uint64_t emission_index = 0;
  for (const double level_photons : o.photons) {
    for (auto &s : sources) {
      s.photons = level_photons / bunch;
      validate(s);
    }
    for (int iev = 0; iev < o.events[0]; ++iev) {
      Event event;
      event.event_number = static_cast<std::int64_t>(events.size()) + 1;
      const auto fire = [&](int flasher, double offset) {
        event.emissions.push_back(
            make_emission(o, sources[flasher], emission_index++, bunch, offset, false));
        event.labels.push_back(flasher + 1);
      };
      if (single >= 0) {
        fire(single, delays[single]);
      } else if (o.fire == "alternate") {
        fire(iev % 4, delays[iev % 4]);
      } else {
        for (int i = 0; i < 4; ++i)
          fire(i, delays[i]);
        if (demo)
          for (int i = 0; i < 4; ++i)
            fire(i, delays[4]);
      }
      events.push_back(std::move(event));
    }
  }
  return events;
}

} // namespace detail

// All events of a run. Throws std::invalid_argument for invalid options.
inline std::vector<Event> build_sources(const Options &o) {
  if (o.photons.empty())
    throw std::invalid_argument("--photons needs at least one value");
  for (const double p : o.photons)
    if (!std::isfinite(p) || p < 0.0)
      throw std::invalid_argument("--photons must be finite and >= 0");
  if (!std::isfinite(o.bunch_size))
    throw std::invalid_argument("--bunchsize must be finite");
  const double bunch = std::max(1.0, o.bunch_size);
  return o.type == Type::dish ? detail::build_dish(o, bunch) : detail::build_secondary(o, bunch);
}

inline void write_stream(std::ostream &out, const Options &o) {
  const auto events = build_sources(o);
  RunSetup setup;
  setup.producer = "obdeect sources flasher";
  setup.run_number = o.run_number;
  setup.observation_level_m = observation_level_m(o);
  setup.atmosphere = o.atmosphere;
  if (o.type == Type::dish) {
    setup.telescope_radius_m = o.camera_radius_cm / 100.0;
    setup.frame = o.frame;
  } else {
    setup.telescope_radius_m = secondary_radius_m;
    setup.frame = obdeect::stream::Frame::telescope;
  }
  programs::write_stream(out, setup, events);
}

} // namespace obdeect::sources::programs::flasher
