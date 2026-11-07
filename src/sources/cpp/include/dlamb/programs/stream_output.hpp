#pragma once

// Shared pieces of the LightEmission ports: option syntax of the original programs,
// the fiducial-sphere rule of LightSource::Emit, and the gstrickd stream writer.
// Header-only, standard library only (plus the dlamb and gstrickd headers).

#include "dlamb/density.hpp"
#include "dlamb/source.hpp"
#include "dlamb/streams.hpp"
#include "dlamb/vec3.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/random.hpp"
#include "gstrickd/stream.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace dlamb::programs {

// Stage id of the program-level random choices (pixel selection). Sources use
// dlamb::stage_sources.
inline constexpr std::uint32_t stage_programs = 0x444C0002U;

// LightSource::Emit only restricts emission to the fiducial sphere when the source is
// farther than 1.02 times the sphere radius (IactLightEmission.cc, rls > 1.02 Rmax).
inline constexpr double fiducial_margin = 1.02;

// FWHM to sigma, as in IactLightEmission.cc.
inline constexpr double sigma_per_fwhm = 2.35482004503095;

// One LightSource::Emit call of the original: a source with its own seed. The bunch size
// multiplies the photon weights and time_offset_ns is the `when` argument of Emit.
struct Emission {
  Source source;
  std::uint64_t seed{};
  double bunch_size{1.0};
  double time_offset_ns{};
};

// One IACT event (one telescope). labels: pixel ids fired (pixled only).
struct Event {
  std::uint32_t number{};
  std::vector<Emission> emissions;
  std::vector<int> labels;
};

// Run and telescope metadata written with every stream.
struct RunSetup {
  std::string producer;
  int run_number{1};
  double observation_level_m{};
  int atmosphere{1};
  std::uint32_t telescope_id{1};
  Vec3 telescope_position{};   // metres, in the photon frame
  double telescope_radius_m{}; // fiducial sphere radius
  gstrickd::Frame frame{gstrickd::Frame::sky};
};

[[noreturn]] inline void bad(const std::string &message) { throw std::invalid_argument(message); }

inline std::string lower(std::string text) {
  for (char &c : text) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return text;
}

inline bool starts_with(const std::string &text, const std::string &prefix) {
  return text.size() >= prefix.size() && text.compare(0, prefix.size(), prefix) == 0;
}

inline std::vector<std::string> split(const std::string &text, char delimiter) {
  std::vector<std::string> parts;
  std::string part;
  for (char c : text) {
    if (c == delimiter) {
      parts.push_back(part);
      part.clear();
    } else {
      part.push_back(c);
    }
  }
  parts.push_back(part);
  return parts;
}

inline double parse_real(const std::string &text, const std::string &what) {
  char *end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  if (text.empty() || end == text.c_str() || *end != '\0' || !std::isfinite(value)) {
    bad(what + ": not a number: '" + text + "'");
  }
  return value;
}

inline long parse_integer(const std::string &text, const std::string &what) {
  char *end = nullptr;
  const long value = std::strtol(text.c_str(), &end, 10);
  if (text.empty() || end == text.c_str() || *end != '\0') {
    bad(what + ": not an integer: '" + text + "'");
  }
  return value;
}

inline std::vector<double> parse_real_list(const std::string &text, const std::string &what) {
  std::vector<double> values;
  for (const auto &part : split(text, ',')) {
    values.push_back(parse_real(part, what));
  }
  return values;
}

inline std::vector<int> parse_integer_list(const std::string &text, const std::string &what) {
  std::vector<int> values;
  for (const auto &part : split(text, ',')) {
    values.push_back(static_cast<int>(parse_integer(part, what)));
  }
  return values;
}

// Value of the option at argv[i] (the next argument); advances i.
inline std::string option_value(int argc, char **argv, int &i) {
  if (i + 1 >= argc) {
    bad(std::string("option ") + argv[i] + " needs a value");
  }
  return argv[++i];
}

// Spectrum as in LightSource::Emit: a wavelength in nm (200-2000) or Cherenkov:from:to
// (1/lambda^2 between from and to; equal limits give a fixed wavelength). Tables read
// from files are not supported.
inline Spectrum parse_spectrum(const std::string &text) {
  if (starts_with(text, "Cherenkov:")) {
    const auto parts = split(text.substr(10), ':');
    if (parts.size() != 2) {
      bad("spectrum '" + text + "': expected Cherenkov:from:to");
    }
    const double from = parse_real(parts[0], "spectrum");
    const double to = parse_real(parts[1], "spectrum");
    if (from > to || from < 200.0 || to > 2000.0) {
      bad("spectrum '" + text + "': need 200 <= from <= to <= 2000 nm");
    }
    if (from == to) {
      return Fixed{from};
    }
    return InverseSquare{from, to};
  }
  // Any plain number (400, 400.5, 4e2) is a wavelength; other text is rejected below.
  char *end = nullptr;
  std::strtod(text.c_str(), &end);
  if (!text.empty() && end != text.c_str() && *end == '\0') {
    const double wavelength = parse_real(text, "spectrum");
    if (wavelength < 200.0 || wavelength > 2000.0) {
      bad("spectrum '" + text + "': wavelength must be 200-2000 nm");
    }
    return Fixed{wavelength};
  }
  bad("spectrum '" + text +
      "': only a wavelength in nm or Cherenkov:from:to is supported (no table files)");
}

// Pulse as in LightSource::Emit: Simple:w[,delay] (or a bare width), Gauss:x (sigma x),
// Gauss-rms:x (sigma x), Gauss-fwhm:x (FWHM x), Exp:tau[,delay], Complex:a,b,c[,delay].
// Widths <= 0 give a time delay only.
inline Pulse parse_pulse(const std::string &text) {
  const std::string key = lower(text);
  if (starts_with(key, "complex:")) {
    const auto values = split(text.substr(8), ',');
    if (values.size() > 4 || values[0].empty()) {
      bad("pulse '" + text + "': expected Complex:flat,gauss,decay[,delay]");
    }
    double v[4] = {0.0, 0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < values.size(); ++i) {
      v[i] = parse_real(values[i], "pulse");
    }
    if (v[0] < 0.0 || v[1] < 0.0 || v[2] < 0.0) {
      bad("pulse '" + text + "': widths must be >= 0");
    }
    return ComplexPulse{v[0], v[1], v[2], v[3]};
  }
  std::string kind = "simple";
  std::string rest = text;
  if (starts_with(key, "simple:")) {
    rest = text.substr(7);
  } else if (starts_with(key, "gauss-fwhm:")) {
    kind = "gauss-fwhm";
    rest = text.substr(11);
  } else if (starts_with(key, "gauss-rms:")) {
    kind = "gauss-rms";
    rest = text.substr(10);
  } else if (starts_with(key, "gauss:")) {
    kind = "gauss-rms";
    rest = text.substr(6);
  } else if (starts_with(key, "exp:")) {
    kind = "exp";
    rest = text.substr(4);
  } else if (text.empty() || !std::isdigit(static_cast<unsigned char>(text[0]))) {
    bad("pulse '" + text +
        "': only Simple:, Gauss:, Gauss-rms:, Gauss-fwhm:, Exp:, Complex: "
        "or a width are supported (no table files)");
  }
  double delay = 0.0;
  const auto comma = rest.find(',');
  if (comma != std::string::npos) {
    delay = parse_real(rest.substr(comma + 1), "pulse delay");
    rest = rest.substr(0, comma);
  }
  const double width = parse_real(rest, "pulse");
  if (kind == "simple") {
    if (width <= 0.0) {
      return Fixed{delay};
    }
    return UniformRange{delay - 0.5 * width, delay + 0.5 * width};
  }
  if (kind == "exp") {
    if (width <= 0.0) {
      return Fixed{delay};
    }
    return Exponential{delay, width};
  }
  // Gauss:x and Gauss-rms:x are sigma x; Gauss-fwhm:x is FWHM x.
  const double sigma = kind == "gauss-fwhm" ? width / sigma_per_fwhm : width;
  if (sigma <= 0.0) {
    return Fixed{delay};
  }
  return GaussianShape{delay, sigma};
}

// Angular profile as in LightSource::Emit: isotropic, Parallel, Cone:deg, FilledCone:deg,
// Gauss:deg / Gauss-rms:deg (2-D rms), Gauss-fwhm:deg, Rayleigh. Tables are not supported.
inline Angular parse_angular(const std::string &text) {
  const std::string key = lower(text);
  constexpr double deg = 3.141592653589793238462643383279502884 / 180.0;
  if (key == "isotropic") {
    return Isotropic{};
  }
  if (key == "parallel") {
    return FixedAngle{0.0};
  }
  if (starts_with(key, "cone:")) {
    return FixedAngle{1.0 - std::cos(deg * parse_real(text.substr(5), "angular"))};
  }
  if (starts_with(key, "filledcone:")) {
    return FilledCone{1.0 - std::cos(deg * parse_real(text.substr(11), "angular"))};
  }
  if (starts_with(key, "gauss-fwhm:")) {
    const double sigma = parse_real(text.substr(11), "angular") / sigma_per_fwhm;
    return TabulatedAngular{gaussian_angular_table(deg * sigma)};
  }
  if (starts_with(key, "gauss-rms:")) {
    return TabulatedAngular{gaussian_angular_table(deg * parse_real(text.substr(10), "angular"))};
  }
  if (starts_with(key, "gauss:")) {
    return TabulatedAngular{gaussian_angular_table(deg * parse_real(text.substr(6), "angular"))};
  }
  if (key == "rayleigh") {
    return TabulatedAngular{rayleigh_angular_table()};
  }
  bad("angular distribution '" + text +
      "': only isotropic, Parallel, Cone:, FilledCone:, Gauss:, Gauss-rms:, Gauss-fwhm:, "
      "Rayleigh are supported (no table files)");
}

// Wavelength range [nm] covered by a spectrum.
inline std::pair<double, double> spectrum_range(const Spectrum &spectrum) {
  return std::visit(
      [](const auto &s) -> std::pair<double, double> {
        using S = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<S, Fixed>) {
          return {s.value, s.value};
        } else if constexpr (std::is_same_v<S, UniformRange> || std::is_same_v<S, InverseSquare>) {
          return {s.low, s.high};
        } else {
          return {s.density.lower(), s.density.upper()};
        }
      },
      spectrum);
}

// The LightSource::Emit rule for the fiducial sphere (radius_m around the origin, the
// telescope centre). Beyond 1.02 R an isotropic emitter is restricted to the cone that
// contains the sphere (SphereTarget, exact); fixed-angle emitters keep all photons (as
// the original, frac = 1); other profiles use a disc target on the sphere's plane (a
// weighted approximation of the original's angular cut).
inline Target fiducial_target(Vec3 position, double radius_m, const Angular &angular) {
  if (!(norm(position) > fiducial_margin * radius_m)) {
    return std::monostate{};
  }
  if (std::holds_alternative<Isotropic>(angular)) {
    // The cone of the original has sin(theta) = R / d (checked against the ff-1m output
    // in scripts/compare_ff-1m.sh); the 1.02 margin only decides when the cut applies.
    return SphereTarget{Vec3{0.0, 0.0, 0.0}, radius_m};
  }
  if (std::holds_alternative<FixedAngle>(angular)) {
    return std::monostate{};
  }
  if (position.z == 0.0) {
    bad("a non-isotropic source in the camera plane outside the fiducial sphere is not supported");
  }
  return DiscTarget{Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, radius_m};
}

// Seed of emission `index` of a program run (splitmix64), so that emissions never share
// random draws.
inline std::uint64_t derive_seed(std::uint64_t base, std::uint64_t index) {
  std::uint64_t z = base + 0x9E3779B97F4A7C15ULL * (index + 1);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

// Uniform deviate of a program-level choice (not a photon draw).
inline double program_uniform(std::uint64_t seed, std::uint32_t event, std::uint64_t draw) {
  return gstrickd::uniform4({seed, stage_programs, event, 0U, draw, 1U})[0];
}

// Writes the complete stream: stream, run, then per event: event, telescope, photon
// chunks, telescope_end, event_end; run_end; stream_end. Photons of one telescope and
// event are numbered consecutively over all emissions (first_index). Sky photons also
// get the zem column (emission altitude a.s.l.).
inline void write_stream(std::ostream &out, const RunSetup &run, const std::vector<Event> &events) {
  gstrickd::Writer writer(out, run.producer);

  double low = std::numeric_limits<double>::infinity();
  double high = -std::numeric_limits<double>::infinity();
  for (const auto &event : events) {
    for (const auto &emission : event.emissions) {
      const auto [a, b] = spectrum_range(emission.source.spectrum);
      low = std::min(low, a);
      high = std::max(high, b);
    }
  }
  const auto frame_value = static_cast<std::int64_t>(run.frame);
  gstrickd::Metadata meta;
  meta["run_number"] = static_cast<std::int64_t>(run.run_number);
  meta["observation_level_m"] = run.observation_level_m;
  meta["n_telescopes"] = std::int64_t{1};
  meta["telescope_x_m"] = std::vector<double>{run.telescope_position.x};
  meta["telescope_y_m"] = std::vector<double>{run.telescope_position.y};
  meta["telescope_z_m"] = std::vector<double>{run.telescope_position.z};
  meta["telescope_r_m"] = std::vector<double>{run.telescope_radius_m};
  meta["atmosphere_id"] = static_cast<std::int64_t>(run.atmosphere);
  meta["refraction_applied"] = std::int64_t{0};
  if (low <= high) {
    meta["wavelength_min_nm"] = low;
    meta["wavelength_max_nm"] = high;
  }
  writer.write_metadata(gstrickd::RecordType::run, meta);

  gstrickd::Metadata telescope;
  telescope["telescope_id"] = static_cast<std::int64_t>(run.telescope_id);
  telescope["frame"] = frame_value;
  telescope["x_m"] = run.telescope_position.x;
  telescope["y_m"] = run.telescope_position.y;
  telescope["z_m"] = run.telescope_position.z;
  telescope["r_m"] = run.telescope_radius_m;

  for (const auto &event : events) {
    gstrickd::Metadata event_meta;
    event_meta["event_number"] = static_cast<std::int64_t>(event.number);
    writer.write_metadata(gstrickd::RecordType::event, event_meta);
    writer.write_metadata(gstrickd::RecordType::telescope, telescope);

    double total_weight = 0.0;
    std::uint64_t index = 0;
    for (const auto &emission : event.emissions) {
      const Source &source = emission.source;
      if (source.frame != run.frame) {
        bad("all emissions of a stream must use the stream frame");
      }
      const std::uint64_t count =
          photon_count(source, emission.seed, event.number, run.telescope_id);
      for (std::uint64_t first = 0; first < count; first += gstrickd::max_chunk_rows) {
        const auto rows = static_cast<std::size_t>(
            std::min<std::uint64_t>(gstrickd::max_chunk_rows, count - first));
        gstrickd::PhotonChunk chunk;
        emit(source, emission.seed, event.number, run.telescope_id, first, rows, chunk);
        chunk.first_index = index + first;
        auto &w = chunk[gstrickd::Column::w];
        auto &t = chunk[gstrickd::Column::t];
        for (std::size_t i = 0; i < rows; ++i) {
          w[i] = static_cast<float>(w[i] * emission.bunch_size);
          t[i] = static_cast<float>(t[i] + emission.time_offset_ns);
          total_weight += w[i];
        }
        if (source.frame == gstrickd::Frame::sky) {
          auto &zem = chunk[gstrickd::Column::zem];
          const auto &z = chunk[gstrickd::Column::z];
          zem.resize(rows);
          for (std::size_t i = 0; i < rows; ++i) {
            zem[i] = static_cast<float>(run.observation_level_m + z[i]);
          }
        }
        writer.write_chunk(chunk);
      }
      index += count;
    }

    gstrickd::Metadata telescope_end;
    telescope_end["telescope_id"] = static_cast<std::int64_t>(run.telescope_id);
    telescope_end["photons"] = total_weight;
    writer.write_metadata(gstrickd::RecordType::telescope_end, telescope_end);
    writer.write_metadata(gstrickd::RecordType::event_end, gstrickd::Metadata{});
  }
  writer.write_metadata(gstrickd::RecordType::run_end, gstrickd::Metadata{});
  writer.finish();
}

// Calls fn(std::ostream&) on stdout (empty path or "-") or on a binary file.
template <class Fn> void write_output(const std::string &path, Fn &&fn) {
  if (path.empty() || path == "-") {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    fn(std::cout);
    std::cout.flush();
    if (!std::cout) {
      bad("write error on stdout");
    }
    return;
  }
  std::ofstream file(path, std::ios::binary);
  if (!file) {
    bad("cannot open output " + path);
  }
  fn(file);
  file.close();
  if (!file) {
    bad("write error on " + path);
  }
}

} // namespace dlamb::programs
