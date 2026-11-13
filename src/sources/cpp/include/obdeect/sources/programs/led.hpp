#pragma once

// LED in front of selected camera pixels, camera frame (docs/spec/tool-led.md).

#include "obdeect/sources/programs/stream_output.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <istream>
#include <map>
#include <optional>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::sources::programs::led {

struct Pixel {
  int id{};
  double x_cm{};
  double y_cm{};
};

struct Camera {
  std::map<int, Pixel> pixels;
  std::vector<std::string> triggers; // pixel-list text of each trigger group
};

struct Options {
  Camera camera;
  std::string use;
  int trigger_random{0};
  bool trigger_all{false};
  int require{0}; // pixels per event, 0 = all of the list
  bool random{false};
  bool slaves{false};
  int events{1};
  double photons{500.0};
  double bunch_size{1.0};
  double camera_radius_cm{123.0};
  std::string spectrum{"400"};
  std::string lightpulse{"Simple:0"};
  std::string angular{"parallel"};
  std::int64_t run_number{1};
  double altitude_m{2150.0};
  std::int64_t atmosphere{1};
  std::uint64_t seed{1};
};

// One pixel of a list: its id, whether it is required, and its pre-sum slaves.
struct ListItem {
  int id{};
  bool required{};
  std::vector<int> slaves;
};

namespace detail {

inline bool is_blank(char c) { return c == ' ' || c == '\t'; }

} // namespace detail

// Pixel-list text of a trigger line: the text after the first "of", with [..] groups removed
// unless slaves are kept.
inline std::string clean_trigger_list(const std::string &line, bool keep_slaves) {
  std::string text = line.substr(0, line.find('\n'));
  const std::size_t of = text.find("of");
  if (of != std::string::npos) {
    std::size_t start = of + 2;
    while (start < text.size() && detail::is_blank(text[start]))
      ++start;
    text = text.substr(start);
  }
  if (keep_slaves)
    return text;
  std::string out;
  bool inside = false;
  for (const char c : text) {
    if (c == '[')
      inside = true;
    else if (c == ']' && inside)
      inside = false;
    else if (!inside)
      out.push_back(c);
  }
  return out;
}

// Reads the pixel and trigger definitions of a camera file.
inline Camera read_camera(std::istream &in, bool keep_slaves) {
  Camera camera;
  std::string line;
  while (std::getline(in, line)) {
    std::size_t start = 0;
    while (start < line.size() && detail::is_blank(line[start]))
      ++start;
    line = line.substr(start);
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream words(line);
    std::string first;
    words >> first;
    const std::string keyword = cli::lower_case(first);
    if (keyword == "pixel") {
      int id = 0, type = 0;
      double x = 0.0, y = 0.0;
      if (!(words >> id >> type >> x >> y)) {
        std::cerr << "led: problem with pixel definition: " << line << '\n';
        continue;
      }
      if (id < 0)
        throw std::invalid_argument("negative pixel id " + std::to_string(id));
      if (!camera.pixels.emplace(id, Pixel{id, x, y}).second)
        throw std::invalid_argument("duplicate pixel id " + std::to_string(id));
    } else if (keyword == "rotate" || keyword == "pixtype") {
      continue;
    } else if (line.find("Trigger") != std::string::npos) {
      camera.triggers.push_back(clean_trigger_list(line, keep_slaves));
    }
  }
  return camera;
}

// Parses a pixel list: [+]id [[a,b,..]] separated by blanks; '%' or a newline ends it.
inline std::vector<ListItem> parse_pixel_list(const std::string &text) {
  std::vector<ListItem> items;
  std::set<int> seen;
  std::size_t i = 0;
  const auto skip_blanks = [&] {
    while (i < text.size() && detail::is_blank(text[i]))
      ++i;
  };
  const auto read_id = [&](int &id) {
    if (i >= text.size() || !std::isdigit(static_cast<unsigned char>(text[i])))
      return false;
    long long v = 0;
    while (i < text.size() && std::isdigit(static_cast<unsigned char>(text[i]))) {
      v = v * 10 + (text[i] - '0');
      if (v > 1000000000LL)
        return false;
      ++i;
    }
    id = static_cast<int>(v);
    return true;
  };
  while (true) {
    skip_blanks();
    if (i >= text.size() || text[i] == '%' || text[i] == '\n')
      break;
    ListItem item;
    if (text[i] == '+') {
      item.required = true;
      ++i;
    }
    if (!read_id(item.id))
      throw std::invalid_argument("invalid in list of pixels to use: '" + text.substr(i) + "'");
    skip_blanks();
    if (i < text.size() && text[i] == '[') {
      ++i;
      while (true) {
        while (i < text.size() && (detail::is_blank(text[i]) || text[i] == ','))
          ++i;
        if (i >= text.size())
          throw std::invalid_argument("pre-sum pixel group not terminated in '" + text + "'");
        if (text[i] == ']') {
          ++i;
          break;
        }
        int slave = 0;
        if (!read_id(slave))
          throw std::invalid_argument("invalid in list of pixels to use: '" + text.substr(i) + "'");
        item.slaves.push_back(slave);
      }
    }
    if (!seen.insert(item.id).second) {
      std::cerr << "duplicated pixel " << item.id << " in use list, used only once\n";
      continue;
    }
    items.push_back(std::move(item));
  }
  return items;
}

namespace detail {

inline void check_options(const Options &o) {
  if (o.camera.pixels.empty())
    throw std::runtime_error("the camera definition has no pixels");
  if (o.events < 0)
    throw std::invalid_argument("--events must be >= 0");
  if (o.require < 0)
    throw std::invalid_argument("--require must be >= 0");
  if (!std::isfinite(o.photons) || o.photons < 0.0)
    throw std::invalid_argument("--photons must be finite and >= 0");
  if (!std::isfinite(o.bunch_size) || !(o.bunch_size > 0.0))
    throw std::invalid_argument("--bunchsize must be finite and > 0");
  if (!std::isfinite(o.camera_radius_cm) || !(o.camera_radius_cm > 0.0))
    throw std::invalid_argument("--camera-radius must be finite and > 0");
  const int modes =
      (!o.use.empty() ? 1 : 0) + (o.trigger_random > 0 ? 1 : 0) + (o.trigger_all ? 1 : 0);
  if (modes != 1)
    throw std::invalid_argument("give exactly one of --use, --use-trg-rand, --use-trg-all");
}

} // namespace detail

// All events of a run. Throws std::invalid_argument for invalid options.
inline std::vector<Event> build_sources(const Options &o) {
  detail::check_options(o);
  const double bunch = std::max(1.0, o.bunch_size);
  const double radius = o.camera_radius_cm / 100.0;
  const Angular angular = parse_angular(o.angular);
  Source base;
  base.frame = obdeect::stream::Frame::camera;
  base.spectrum = parse_spectrum(o.spectrum);
  base.pulse = parse_pulse(o.lightpulse);
  base.photons = o.photons / bunch;

  // Pixel lists.
  std::vector<std::vector<ListItem>> lists;
  if (o.trigger_all) {
    if (o.camera.triggers.empty())
      throw std::invalid_argument("the camera definition has no trigger groups");
    for (const auto &text : o.camera.triggers)
      lists.push_back(parse_pixel_list(text));
  } else if (o.trigger_random > 0) {
    if (o.camera.triggers.empty())
      throw std::invalid_argument("the camera definition has no trigger groups");
    const auto groups = o.camera.triggers.size();
    std::uint64_t draw_index = 0;
    for (int n = 0; n < o.trigger_random; ++n) {
      std::size_t j;
      do {
        j = static_cast<std::size_t>(
            std::floor(static_cast<double>(groups) * program_uniform(o.seed, 0, draw_index++)));
      } while (j >= groups);
      lists.push_back(parse_pixel_list(o.camera.triggers[j]));
    }
  } else {
    lists.push_back(parse_pixel_list(clean_trigger_list(o.use, o.slaves)));
  }

  std::vector<Event> events;
  std::uint64_t emission_index = 0;
  const auto add_event = [&](const std::vector<const ListItem *> &lit) {
    Event event;
    event.event_number = static_cast<std::int64_t>(events.size()) + 1;
    const auto emit_pixel = [&](int id) {
      const Pixel &pixel = o.camera.pixels.at(id);
      Source s = base;
      const Vec3 position{pixel.x_cm / 100.0, pixel.y_cm / 100.0, 0.0};
      s.geometry = Emitter{position, {0.0, 0.0, -1.0}, angular};
      s.target = fiducial_target(position, radius, angular);
      validate(s);
      Emission e;
      e.source = s;
      e.seed = derive_seed(o.seed, emission_index++);
      e.telescope = 1;
      e.bunch_size = bunch;
      event.emissions.push_back(std::move(e));
    };
    for (const ListItem *item : lit) {
      emit_pixel(item->id);
      event.labels.push_back(item->id);
      if (o.slaves)
        for (const int slave : item->slaves)
          emit_pixel(slave);
    }
    events.push_back(std::move(event));
  };

  for (const auto &list : lists) {
    if (list.empty())
      throw std::invalid_argument("failed to obtain a pixel list");
    for (const auto &item : list) {
      if (!o.camera.pixels.count(item.id))
        throw std::invalid_argument("pixel ID " + std::to_string(item.id) +
                                    " is not defined in the camera file");
      if (o.slaves)
        for (const int slave : item.slaves)
          if (!o.camera.pixels.count(slave))
            throw std::invalid_argument("pixel ID " + std::to_string(slave) +
                                        " is not defined in the camera file");
    }
    const std::size_t size = list.size();
    const std::size_t nreq = o.require == 0 ? size : static_cast<std::size_t>(o.require);
    if (nreq > size) {
      std::cerr << "led: required " << nreq << " pixels but the list has " << size
                << "; list skipped\n";
      continue;
    }
    std::vector<const ListItem *> all;
    for (const auto &item : list)
      all.push_back(&item);
    for (int repeat = 0; repeat < o.events; ++repeat) {
      if (nreq == size) {
        add_event(all);
      } else if (o.random) {
        const auto event_number = static_cast<std::uint32_t>(events.size() + 1);
        std::vector<bool> used(size, false);
        std::vector<const ListItem *> order; // fired in draw order
        std::uint64_t draw_index = 0;
        while (order.size() < nreq) {
          const auto p = static_cast<std::size_t>(std::floor(
              static_cast<double>(size) * program_uniform(o.seed, event_number, draw_index++)));
          if (p >= size || used[p])
            continue;
          used[p] = true;
          order.push_back(&list[p]);
        }
        add_event(order);
      } else {
        std::vector<std::size_t> required, optional;
        for (std::size_t k = 0; k < size; ++k)
          (list[k].required ? required : optional).push_back(k);
        const auto emit_combination = [&](const std::vector<std::size_t> &positions) {
          std::vector<std::size_t> sorted = positions;
          std::sort(sorted.begin(), sorted.end());
          std::vector<const ListItem *> lit;
          for (const std::size_t k : sorted)
            lit.push_back(&list[k]);
          add_event(lit);
        };
        if (required.size() >= nreq) {
          emit_combination(required);
        } else {
          const std::size_t extra = nreq - required.size();
          std::vector<std::size_t> pick(extra);
          for (std::size_t k = 0; k < extra; ++k)
            pick[k] = k;
          while (true) {
            std::vector<std::size_t> positions = required;
            for (const std::size_t k : pick)
              positions.push_back(optional[k]);
            emit_combination(positions);
            // Next combination of `extra` of optional.size() in lexicographic order.
            std::size_t k = extra;
            while (k > 0 && pick[k - 1] == optional.size() - extra + k - 1)
              --k;
            if (k == 0)
              break;
            ++pick[k - 1];
            for (std::size_t m = k; m < extra; ++m)
              pick[m] = pick[m - 1] + 1;
          }
        }
      }
    }
  }
  return events;
}

inline void write_stream(std::ostream &out, const Options &o) {
  const auto events = build_sources(o);
  RunSetup setup;
  setup.producer = "obdeect sources led";
  setup.run_number = o.run_number;
  setup.observation_level_m = o.altitude_m;
  setup.atmosphere = o.atmosphere;
  setup.telescope_radius_m = o.camera_radius_cm / 100.0;
  setup.frame = obdeect::stream::Frame::camera;
  programs::write_stream(out, setup, events);
}

// The `--list-only` text: one line per event with the sorted pixel labels.
inline std::string label_lines(const std::vector<Event> &events) {
  std::string text;
  for (const auto &event : events) {
    std::vector<int> labels = event.labels;
    std::sort(labels.begin(), labels.end());
    text += "(@) ";
    for (const int label : labels)
      text += " " + std::to_string(label);
    text += "\n";
  }
  return text;
}

} // namespace obdeect::sources::programs::led
