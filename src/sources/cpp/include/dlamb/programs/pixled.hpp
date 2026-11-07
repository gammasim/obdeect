#pragma once

// Port of sim_telarray/LightEmission/pixled.cc: an LED in front of selected camera pixels,
// to test trigger definitions.
//
// Geometry: LED at the pixel centre (x, y) of the camera definition, z = 0 (camera plane,
// dist_cm = 0 in pixled.cc:508), pointing down (parallel to -z, pixled.cc:511, 757). The
// photons are already on the camera plane, so they are written in the camera frame C,
// in the pixel-file coordinates (no y flip: sim_telarray and ps2simtel apply it).

#include "dlamb/programs/stream_output.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <istream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace dlamb::programs::pixled {

// Pixel centre in the camera frame [cm].
struct PixelPosition {
  double x_cm{};
  double y_cm{};
};

struct Camera {
  // Keyed by pixel ID as cam[id] (pixled.cc:175); the file order does not matter.
  std::map<int, PixelPosition> pixels;
  std::vector<std::string> triggers; // cleaned trigger-group pixel lists
};

// One pixel of a pixel list: '+' marks a required pixel, [..] lists pre-sum slaves.
struct Use {
  int id{};
  bool required{};
  std::vector<int> slaves;
};

// Options with the names of pixled.cc. Defaults are the values in the code (the --help
// text says 2500000 photons; the code uses 500, pixled.cc:509).
struct Options {
  Camera camera;                      // read from --camera-file (pixled.cc:652-667)
  bool slaves{false};                 // --slaves, pixled.cc:540-545
  std::string use;                    // --use pixel-list, pixled.cc:583-590
  int trigger_random{0};              // --use-trg-rand n, pixled.cc:599-607 (0 = off)
  bool trigger_all{false};            // --use-trg-all, pixled.cc:608-617
  int require{0};                     // --require n, pixled.cc:618-619 (0 = all pixels)
  bool random{false};                 // --random n, pixled.cc:620-624 (uses --require n)
  int events{1};                      // --events, pixled.cc:629-630
  double photons{500.0};              // --photons [per pixel], pixled.cc:509, 571-572
  double bunch_size{1.0};             // --bunchsize, pixled.cc:510, 573-574
  double camera_radius_cm{123.0};     // --camera-radius [cm], pixled.cc:507, 575-576
  std::string spectrum{"400"};        // --spectrum, pixled.cc:511, 579-580
  std::string lightpulse{"Simple:0"}; // --lightpulse, pixled.cc:511, 581-582
  std::string angular{"parallel"};    // pixled.cc:511 (fixed to Parallel in the original)
  int run_number{1};                  // --run, pixled.cc:512, 625-626
  double altitude_m{2150.0};          // --altitude [m], pixled.cc:506, 569-570
  int atmosphere{1};                  // --atmosphere, pixled.cc:512, 627-628 (numbers only)
  std::uint64_t seed{1};              // new: Philox seed (the original uses the global RNG)
};

// Reads the Pixel and Trigger lines of a sim_telarray camera definition (pixled.cc:139-199).
inline Camera read_camera(std::istream &in, bool keep_slaves);

// Trigger or pixel list without the "...of N:" prefix and, unless keep_slaves, without
// the [..] slave lists (pixled.cc:96-131). The original's copy of the character after a
// removed [..] group is dropped here; no list depends on it.
inline std::string clean_trigger_list(const std::string &text, bool keep_slaves) {
  const std::string line = text.substr(0, text.find('\n'));
  std::size_t start = 0;
  const auto of = line.find("of");
  if (of != std::string::npos) {
    start = of + 2;
    while (start < line.size() && (line[start] == ' ' || line[start] == '\t')) {
      ++start;
    }
  }
  std::string out;
  for (std::size_t i = start; i < line.size(); ++i) {
    if (line[i] == '[' && !keep_slaves) {
      while (i < line.size() && line[i] != ']') {
        ++i;
      }
      continue;
    }
    out.push_back(line[i]);
  }
  return out;
}

inline Camera read_camera(std::istream &in, bool keep_slaves) {
  Camera camera;
  std::string line;
  while (std::getline(in, line)) {
    const auto first = line.find_first_not_of(" \t");
    if (first == std::string::npos || line[first] == '#') {
      continue;
    }
    const std::string body = line.substr(first);
    const std::string key = lower(body.substr(0, 8));
    if (starts_with(key, "pixel") && body.size() > 5 && (body[5] == ' ' || body[5] == '\t')) {
      std::istringstream fields(body.substr(5));
      int id = 0;
      int type = 0;
      double x = 0.0;
      double y = 0.0;
      if (!(fields >> id >> type >> x >> y)) {
        std::cerr << "pixled: problem with pixel definition: " << line << "\n";
        continue;
      }
      if (id < 0) {
        bad("pixled: negative pixel ID in camera definition: " + line);
      }
      if (!camera.pixels.emplace(id, PixelPosition{x, y}).second) {
        bad("pixled: duplicate pixel ID " + std::to_string(id) + " in camera definition");
      }
    } else if (starts_with(key, "rotate") || starts_with(key, "pixtype")) {
      continue;
    } else if (line.find("Trigger") != std::string::npos) {
      camera.triggers.push_back(clean_trigger_list(line, keep_slaves));
    }
  }
  return camera;
}

namespace detail {

inline int read_int(const std::string &text, std::size_t &i) {
  int value = 0;
  while (i < text.size() && std::isdigit(static_cast<unsigned char>(text[i]))) {
    value = value * 10 + (text[i] - '0');
    ++i;
  }
  return value;
}

inline bool is_digit_at(const std::string &text, std::size_t i) {
  return i < text.size() && std::isdigit(static_cast<unsigned char>(text[i])) != 0;
}

} // namespace detail

// Pixel list such as "3 +5 7[8,9]": pixel ids separated by blanks, '+' for required
// pixels, [..] for pre-sum slaves, '%' starts a comment (pixled.cc:205-334, simplified:
// no separate fixed-width buffers). Duplicate ids are skipped with a warning.
inline std::vector<Use> parse_pixel_list(const std::string &text) {
  std::vector<Use> result;
  std::size_t i = 0;
  const auto skip_blanks = [&] {
    while (i < text.size() && (text[i] == ' ' || text[i] == '\t')) {
      ++i;
    }
  };
  for (;;) {
    skip_blanks();
    if (i >= text.size() || text[i] == '%' || text[i] == '\n') {
      break;
    }
    Use use;
    if (text[i] == '+') {
      use.required = true;
      ++i;
    }
    if (!detail::is_digit_at(text, i)) {
      bad("pixled: invalid in list of pixels to use: '" + text.substr(i) + "'");
    }
    use.id = detail::read_int(text, i);
    skip_blanks();
    if (i < text.size() && text[i] == '[') {
      ++i;
      for (;;) {
        skip_blanks();
        while (i < text.size() && text[i] == ',') {
          ++i;
          skip_blanks();
        }
        if (i < text.size() && text[i] == ']') {
          ++i;
          break;
        }
        if (!detail::is_digit_at(text, i)) {
          bad("pixled: pre-sum pixel group not terminated in '" + text + "'");
        }
        use.slaves.push_back(detail::read_int(text, i));
      }
    }
    const bool duplicate =
        std::any_of(result.begin(), result.end(), [&](const Use &u) { return u.id == use.id; });
    if (duplicate) {
      std::cerr << "pixled: duplicated pixel " << use.id << " in use list, used only once\n";
      continue;
    }
    result.push_back(use);
  }
  return result;
}

inline void validate(const Options &o) {
  const int modes =
      (o.use.empty() ? 0 : 1) + (o.trigger_random > 0 ? 1 : 0) + (o.trigger_all ? 1 : 0);
  if (modes != 1) {
    bad("pixled: give exactly one of --use, --use-trg-rand n, --use-trg-all");
  }
  if (o.camera.pixels.empty()) {
    bad("pixled: no pixels available in the camera definition");
  }
  if (o.events < 0 || o.require < 0) {
    bad("pixled: --events and --require must be >= 0");
  }
  if (!(o.photons >= 0.0) || !(o.bunch_size > 0.0) || !(o.camera_radius_cm > 0.0)) {
    bad("pixled: photons >= 0, bunch size > 0 and camera radius > 0 required");
  }
}

// Recursive selection of pixel combinations (pixled.cc:424-497): the required pixels plus
// the next combination of available ones, in increasing order.
template <class Callback>
void for_each_combination(const std::vector<std::size_t> &selected,
                          const std::vector<std::size_t> &available, std::size_t nreq,
                          Callback &callback) {
  if (selected.size() >= nreq) {
    callback(selected);
    return;
  }
  for (std::size_t k = 0; k < available.size(); ++k) {
    auto next = selected;
    next.push_back(available[k]);
    const std::vector<std::size_t> rest(available.begin() + static_cast<std::ptrdiff_t>(k) + 1,
                                        available.end());
    for_each_combination(next, rest, nreq, callback);
  }
}

// Sources of all events (pixled.cc:679-916). Pixel lists come from --use, or from trigger
// groups of the camera file (--use-trg-all, --use-trg-rand n). Per list and per event:
// all pixels (--require 0), or --random n distinct random pixels, or every combination of
// --require n pixels as one event each. Each lit pixel is one emission from its centre.
inline std::vector<Event> build_sources(const Options &o) {
  validate(o);
  const double bunch = o.bunch_size < 1.0 ? 1.0 : o.bunch_size; // IactLightEmission.cc:1514
  const double radius_m = 0.01 * o.camera_radius_cm;
  const Angular angular = parse_angular(o.angular);
  const Spectrum spectrum = parse_spectrum(o.spectrum);
  const Pulse pulse = parse_pulse(o.lightpulse);

  std::vector<std::vector<Use>> lists;
  if (o.trigger_all) {
    if (o.camera.triggers.empty()) {
      bad("pixled: the camera definition has no trigger groups");
    }
    for (const auto &group : o.camera.triggers) {
      lists.push_back(parse_pixel_list(group));
    }
  } else if (o.trigger_random > 0) {
    if (o.camera.triggers.empty()) {
      bad("pixled: the camera definition has no trigger groups");
    }
    const auto groups = o.camera.triggers.size();
    std::uint64_t draw = 0;
    for (int k = 0; k < o.trigger_random; ++k) {
      std::size_t j = groups;
      while (j >= groups) {
        j = static_cast<std::size_t>(
            std::floor(static_cast<double>(groups) * program_uniform(o.seed, 0, draw++)));
      }
      lists.push_back(parse_pixel_list(o.camera.triggers[j]));
    }
  } else {
    lists.push_back(parse_pixel_list(clean_trigger_list(o.use, o.slaves)));
  }

  std::vector<Event> result;
  std::uint64_t emission = 0;
  const auto add_pixel = [&](Event &event, int id, bool label) {
    const auto found = o.camera.pixels.find(id);
    if (found == o.camera.pixels.end()) {
      bad("pixled: pixel ID " + std::to_string(id) + " is not defined in the camera file");
    }
    const auto &pixel = found->second;
    const Vec3 position{0.01 * pixel.x_cm, 0.01 * pixel.y_cm, 0.0};
    Source source;
    source.frame = gstrickd::Frame::camera;
    source.geometry = Emitter{position, Vec3{0.0, 0.0, -1.0}, angular};
    source.spectrum = spectrum;
    source.pulse = pulse;
    source.photons = o.photons / bunch;
    source.target = fiducial_target(position, radius_m, angular);
    event.emissions.push_back({source, derive_seed(o.seed, emission++), bunch, 0.0});
    if (label) {
      event.labels.push_back(id);
    }
  };
  const auto add_use = [&](Event &event, const Use &use) {
    add_pixel(event, use.id, true);
    if (o.slaves) {
      for (const int slave : use.slaves) {
        add_pixel(event, slave, false);
      }
    }
  };

  for (const auto &list : lists) {
    if (list.empty()) {
      bad("pixled: failed to obtain a pixel list");
    }
    std::size_t nreq = o.require == 0 ? list.size() : static_cast<std::size_t>(o.require);
    if (nreq > list.size()) {
      std::cerr << "pixled: required " << nreq << " pixels but the list has " << list.size()
                << "; list skipped\n";
      continue;
    }
    for (int iev = 0; iev < o.events; ++iev) {
      Event event;
      event.number = static_cast<std::uint32_t>(result.size() + 1);
      if (nreq == list.size()) {
        for (const auto &use : list) {
          add_use(event, use);
        }
        result.push_back(std::move(event));
      } else if (o.random) {
        std::uint64_t draw = 0;
        std::vector<int> fired;
        while (fired.size() < nreq) {
          const auto p = static_cast<std::size_t>(std::floor(
              static_cast<double>(list.size()) * program_uniform(o.seed, event.number, draw++)));
          if (p >= list.size()) {
            continue;
          }
          const bool used = std::find(fired.begin(), fired.end(), list[p].id) != fired.end();
          if (used) {
            continue;
          }
          fired.push_back(list[p].id);
          add_use(event, list[p]);
        }
        result.push_back(std::move(event));
      } else {
        std::vector<std::size_t> selected;
        std::vector<std::size_t> available;
        for (std::size_t j = 0; j < list.size(); ++j) {
          (list[j].required ? selected : available).push_back(j);
        }
        const auto emit_combination = [&](const std::vector<std::size_t> &combination) {
          Event combo;
          combo.number = static_cast<std::uint32_t>(result.size() + 1);
          for (const std::size_t p : combination) {
            add_use(combo, list[p]);
          }
          result.push_back(std::move(combo));
        };
        for_each_combination(selected, available, nreq, emit_combination);
      }
    }
  }
  return result;
}

inline RunSetup run_setup(const Options &o) {
  RunSetup run;
  run.producer = "dlamb-pixled";
  run.run_number = o.run_number;
  run.observation_level_m = o.altitude_m;
  run.atmosphere = o.atmosphere;
  run.telescope_radius_m = 0.01 * o.camera_radius_cm;
  run.frame = gstrickd::Frame::camera;
  return run;
}

// Pixel labels of an event in the --list-only format of pixled.cc:443-453.
inline void print_labels(std::ostream &out, const Event &event) {
  auto sorted = event.labels;
  std::sort(sorted.begin(), sorted.end());
  out << "(@) ";
  for (const int id : sorted) {
    out << ' ' << id;
  }
  out << '\n';
}

inline void write_stream(std::ostream &out, const Options &o) {
  dlamb::programs::write_stream(out, run_setup(o), build_sources(o));
}

} // namespace dlamb::programs::pixled
