#pragma once

// Shared parts of the flasher and led programs: option-text parsers, the fiducial-sphere rule,
// and the program stream writer (docs/spec/stream-writer.md, writer A).

#include "obdeect/sources/cli.hpp"
#include "obdeect/sources/density.hpp"
#include "obdeect/sources/emission_stream.hpp"
#include "obdeect/sources/source.hpp"
#include "obdeect/sources/vec3.hpp"

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"
#include "obdeect/stream/stages.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::sources::programs {

// Philox stage of program-level choices (obdeect::stream::stage_programs).
inline constexpr std::uint32_t stage_programs = obdeect::stream::stage_programs;

// Seed of emission i of a program (same function as obdeect::sources::source_seed).
inline std::uint64_t derive_seed(std::uint64_t base, std::uint64_t i) {
  return source_seed(base, i);
}

// Uniform deviate of a program-level choice.
inline double program_uniform(std::uint64_t seed, std::uint32_t event, std::uint64_t draw_index) {
  return obdeect::stream::uniform4({seed, stage_programs, event, 0, draw_index, 1})[0];
}

// One source emitted with a bunch size and a time offset; to_camera_plane moves the photons in a
// straight line to the plane z = 0.
struct Emission {
  Source source;
  std::uint64_t seed{1};
  std::uint32_t telescope{1};
  double bunch_size{1.0};
  double time_offset_ns{};
  bool to_camera_plane{};
};

struct Event {
  std::int64_t event_number{1};
  std::vector<Emission> emissions;
  std::vector<int> labels; // pixels or flashers that fire (informational)
};

struct RunSetup {
  std::string producer;
  std::int64_t run_number{1};
  double observation_level_m{2150.0};
  std::int64_t atmosphere{1};
  std::uint32_t telescope_id{1};
  Vec3 telescope_position{};
  double telescope_radius_m{};
  obdeect::stream::Frame frame{obdeect::stream::Frame::sky};
};

// Parsers of option text (messages are `<option>: ...`).

inline double parse_real(const std::string &text, const std::string &option) {
  return cli::parse_number(text, option);
}

inline long long parse_integer_value(const std::string &text, const std::string &option) {
  try {
    return cli::parse_integer(text, option);
  } catch (const std::invalid_argument &) {
    throw std::invalid_argument(option + ": not an integer: '" + text + "'");
  }
}

inline std::vector<double> parse_real_list(const std::string &text, const std::string &option) {
  std::vector<double> values;
  for (const auto &item : cli::split_items(text, option))
    values.push_back(parse_real(item, option));
  return values;
}

inline Spectrum parse_spectrum(const std::string &text) {
  const auto in_range = [&](double v) {
    if (!(v >= 200.0 && v <= 2000.0))
      throw std::invalid_argument("spectrum '" + text + "': wavelengths must lie in 200..2000 nm");
    return v;
  };
  if (cli::is_number(text))
    return Fixed{in_range(parse_real(text, "spectrum"))};
  if (cli::starts_with_nocase(text, "cherenkov:")) {
    const std::string rest = text.substr(10);
    const std::size_t colon = rest.find(':');
    if (colon == std::string::npos || rest.find(':', colon + 1) != std::string::npos)
      throw std::invalid_argument("spectrum '" + text + "': Cherenkov needs two fields A:B");
    const double a = in_range(parse_real(rest.substr(0, colon), "spectrum"));
    const double b = in_range(parse_real(rest.substr(colon + 1), "spectrum"));
    if (a > b)
      throw std::invalid_argument("spectrum '" + text + "': Cherenkov range needs A <= B");
    if (a == b)
      return Fixed{a};
    return InverseSquare{a, b};
  }
  throw std::invalid_argument("spectrum '" + text +
                              "': supported are a wavelength in nm or Cherenkov:A:B");
}

inline Pulse parse_pulse(const std::string &text) {
  const std::string supported = "supported are Simple:W[,D], W, Gauss:S[,D], Gauss-rms:S[,D], "
                                "Gauss-fwhm:F[,D], Exp:T[,D], Complex:A,B,C[,D]";
  const auto bad = [&](const std::string &why) {
    return std::invalid_argument("pulse '" + text + "': " + why);
  };
  const auto simple = [](double w, double d) {
    return w > 0.0 ? Pulse{UniformRange{d - 0.5 * w, d + 0.5 * w}} : Pulse{Fixed{d}};
  };
  if (!text.empty() && std::isdigit(static_cast<unsigned char>(text[0]))) {
    if (!cli::is_number(text))
      throw bad(supported);
    return simple(parse_real(text, "pulse"), 0.0);
  }
  const std::size_t colon = text.find(':');
  if (colon == std::string::npos)
    throw bad(supported);
  const std::string name = cli::lower_case(text.substr(0, colon));
  const std::string rest = text.substr(colon + 1);
  if (name == "complex") {
    std::vector<double> v(4, 0.0);
    std::vector<std::string> items;
    std::size_t start = 0;
    while (true) {
      const std::size_t comma = rest.find(',', start);
      items.push_back(rest.substr(start, comma == std::string::npos ? comma : comma - start));
      if (comma == std::string::npos)
        break;
      start = comma + 1;
    }
    if (items.size() > 4 || items[0].empty())
      throw bad("Complex takes one to four values, the first not empty");
    for (std::size_t i = 0; i < items.size(); ++i)
      if (!items[i].empty())
        v[i] = parse_real(items[i], "pulse");
    if (v[0] < 0.0 || v[1] < 0.0 || v[2] < 0.0)
      throw bad("Complex widths must be >= 0");
    return ComplexPulse{v[0], v[1], v[2], v[3]};
  }
  if (name != "simple" && name != "gauss" && name != "gauss-rms" && name != "gauss-fwhm" &&
      name != "exp")
    throw bad(supported);
  const auto values = parse_real_list(rest, "pulse");
  if (values.size() > 2)
    throw bad("at most two values (width and delay)");
  const double delay = values.size() > 1 ? values[1] : 0.0;
  const double s = values[0];
  if (name == "simple")
    return simple(s, delay);
  if (name == "exp")
    return s > 0.0 ? Pulse{Exponential{delay, s}} : Pulse{Fixed{delay}};
  const double sigma = name == "gauss-fwhm" ? s / fwhm_per_sigma : s;
  return sigma > 0.0 ? Pulse{GaussianShape{delay, sigma}} : Pulse{Fixed{delay}};
}

inline Angular parse_angular(const std::string &text) {
  const std::string lower = cli::lower_case(text);
  if (lower == "isotropic")
    return Isotropic{};
  if (lower == "parallel")
    return FixedAngle{0.0};
  if (lower == "rayleigh")
    return TabulatedAngular{rayleigh_angular_table()};
  const std::size_t colon = text.find(':');
  if (colon != std::string::npos) {
    const std::string name = cli::lower_case(text.substr(0, colon));
    const double value = parse_real(text.substr(colon + 1), "angular distribution");
    const auto degrees = [&] {
      if (!(value >= 0.0 && value <= 180.0))
        throw std::invalid_argument("angular distribution '" + text +
                                    "': angle must be finite and in [0, 180] deg");
      return value * pi / 180.0;
    };
    if (name == "cone")
      return FixedAngle{1.0 - std::cos(degrees())};
    if (name == "filledcone") {
      if (!(value > 0.0))
        throw std::invalid_argument("angular distribution '" + text + "': angle must be > 0");
      return FilledCone{1.0 - std::cos(degrees())};
    }
    if (name == "gauss" || name == "gauss-rms")
      return TabulatedAngular{gaussian_angular_table(value * pi / 180.0)};
    if (name == "gauss-fwhm")
      return TabulatedAngular{gaussian_angular_table(value / fwhm_per_sigma * pi / 180.0)};
  }
  throw std::invalid_argument("angular distribution '" + text +
                              "': supported are isotropic, Parallel, Cone:A, FilledCone:A, "
                              "Gauss:S, Gauss-fwhm:F, Rayleigh");
}

// Target of a source in front of a sphere of radius R about the origin: none inside 1.02 R, a
// cone target for isotropic emission, none for a fixed angle, otherwise a disc weighted target.
inline Target fiducial_target(const Vec3 &position, double radius_m, const Angular &angular) {
  if (norm(position) <= 1.02 * radius_m)
    return std::monostate{};
  if (std::holds_alternative<Isotropic>(angular))
    return SphereTarget{{0.0, 0.0, 0.0}, radius_m};
  if (std::holds_alternative<FixedAngle>(angular))
    return std::monostate{};
  if (position.z == 0.0)
    throw std::invalid_argument("a non-isotropic source in the camera plane outside the fiducial "
                                "sphere is not supported");
  return DiscTarget{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, radius_m};
}

// Writes the run, then each event with one telescope block holding all its emissions.
inline void write_stream(std::ostream &out, const RunSetup &setup,
                         const std::vector<Event> &events) {
  using obdeect::stream::Column;
  obdeect::stream::Writer writer(out, setup.producer);
  obdeect::stream::Metadata run;
  run["run_number"] = setup.run_number;
  run["observation_level_m"] = setup.observation_level_m;
  run["n_telescopes"] = std::int64_t{1};
  run["telescope_x_m"] = std::vector<double>{setup.telescope_position.x};
  run["telescope_y_m"] = std::vector<double>{setup.telescope_position.y};
  run["telescope_z_m"] = std::vector<double>{setup.telescope_position.z};
  run["telescope_r_m"] = std::vector<double>{setup.telescope_radius_m};
  run["atmosphere_id"] = setup.atmosphere;
  run["refraction_applied"] = std::int64_t{0};
  bool any = false;
  double wl_min = 0.0, wl_max = 0.0;
  for (const auto &event : events)
    for (const auto &emission : event.emissions) {
      const auto [lo, hi] = spectrum_range(emission.source.spectrum);
      wl_min = any ? std::min(wl_min, lo) : lo;
      wl_max = any ? std::max(wl_max, hi) : hi;
      any = true;
    }
  if (any) {
    run["wavelength_min_nm"] = wl_min;
    run["wavelength_max_nm"] = wl_max;
  }
  writer.write_metadata(obdeect::stream::RecordType::run, run);

  obdeect::stream::PhotonChunk chunk;
  for (const auto &event : events) {
    const auto event_key = static_cast<std::uint32_t>(event.event_number);
    writer.write_metadata(obdeect::stream::RecordType::event,
                          {{"event_number", event.event_number}});
    obdeect::stream::Metadata tm;
    tm["telescope_id"] = static_cast<std::int64_t>(setup.telescope_id);
    tm["frame"] = static_cast<std::int64_t>(setup.frame);
    tm["x_m"] = setup.telescope_position.x;
    tm["y_m"] = setup.telescope_position.y;
    tm["z_m"] = setup.telescope_position.z;
    tm["r_m"] = setup.telescope_radius_m;
    writer.write_metadata(obdeect::stream::RecordType::telescope, tm);
    std::uint64_t written = 0;
    double weight = 0.0;
    for (const auto &emission : event.emissions) {
      if (emission.source.frame != setup.frame)
        throw std::invalid_argument("all emissions of a stream must use the stream frame");
      const std::uint64_t count =
          photon_count(emission.source, emission.seed, event_key, setup.telescope_id);
      for (std::uint64_t done = 0; done < count;) {
        const auto rows = static_cast<std::size_t>(
            std::min<std::uint64_t>(obdeect::stream::max_chunk_rows, count - done));
        emit(emission.source, emission.seed, event_key, setup.telescope_id, done, rows, chunk);
        chunk.first_index = written;
        auto &x = chunk[Column::x];
        auto &y = chunk[Column::y];
        auto &z = chunk[Column::z];
        auto &t = chunk[Column::t];
        auto &w = chunk[Column::w];
        if (emission.to_camera_plane) {
          for (std::size_t i = 0; i < rows; ++i) {
            const double cz = chunk[Column::cz][i];
            if (!(cz < 0.0) || !(z[i] >= 0.0F))
              throw std::invalid_argument("a photon does not travel towards the camera plane");
            const double path = static_cast<double>(z[i]) / -cz;
            x[i] = static_cast<float>(x[i] + path * chunk[Column::cx][i]);
            y[i] = static_cast<float>(y[i] + path * chunk[Column::cy][i]);
            t[i] = static_cast<float>(t[i] + path / speed_of_light_m_per_ns);
            z[i] = 0.0F;
          }
        }
        for (std::size_t i = 0; i < rows; ++i) {
          w[i] = static_cast<float>(static_cast<double>(w[i]) * emission.bunch_size);
          t[i] = static_cast<float>(static_cast<double>(t[i]) + emission.time_offset_ns);
          weight += static_cast<double>(w[i]);
        }
        if (setup.frame == obdeect::stream::Frame::sky) {
          auto &zem = chunk[Column::zem];
          zem.resize(rows);
          for (std::size_t i = 0; i < rows; ++i)
            zem[i] = static_cast<float>(setup.observation_level_m + static_cast<double>(z[i]));
        }
        writer.write_chunk(chunk);
        written += rows;
        done += rows;
      }
    }
    writer.write_metadata(
        obdeect::stream::RecordType::telescope_end,
        {{"telescope_id", static_cast<std::int64_t>(setup.telescope_id)}, {"photons", weight}});
    writer.write_metadata(obdeect::stream::RecordType::event_end, {});
  }
  writer.write_metadata(obdeect::stream::RecordType::run_end, {});
  writer.finish();
}

} // namespace obdeect::sources::programs
