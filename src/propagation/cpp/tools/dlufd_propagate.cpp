// dlufd-propagate: the atmosphere stage on a gstrickd photon stream (docs/propagation.md).
// Reads a stream from a file or stdin, writes the processed stream to a file or stdout, and
// prints a summary to stderr.

#include "dlufd/propagator.hpp"

using namespace dlufd;

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

void usage() {
  std::fprintf(stderr,
               "usage: dlufd-propagate --transmission FILE --site-altitude M [--profile FILE]\n"
               "                       [--envelope NAME=FILE]... [--wavelength-range LO HI]\n"
               "                       --seed N [--mode thinned|weighted] [--store-p]\n"
               "                       [--scattering FILE] [--scattering-points N]\n"
               "                       [-i FILE|-] [-o FILE|-]\n");
}

// Thrown for bad command lines (exit status 2).
struct UsageError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

double parse_double(const std::string &text, const char *what) {
  char *end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  if (text.empty() || end != text.c_str() + text.size()) {
    throw UsageError(std::string("invalid number for ") + what + ": " + text);
  }
  return value;
}

std::uint32_t parse_points(const std::string &text) {
  char *end = nullptr;
  const long value = std::strtol(text.c_str(), &end, 10);
  if (text.empty() || end != text.c_str() + text.size() || value < 1 ||
      value > static_cast<long>(max_points_per_photon)) {
    throw UsageError("--scattering-points needs an integer from 1 to 1024: " + text);
  }
  return static_cast<std::uint32_t>(value);
}

std::uint64_t parse_seed(const std::string &text) {
  char *end = nullptr;
  errno = 0;
  const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
  if (text.empty() || text[0] == '-' || end != text.c_str() + text.size() || errno != 0) {
    throw UsageError("invalid seed: " + text);
  }
  return static_cast<std::uint64_t>(value);
}

void set_binary_stdio() {
#ifdef _WIN32
  _setmode(_fileno(stdin), _O_BINARY);
  _setmode(_fileno(stdout), _O_BINARY);
#endif
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  set_binary_stdio();
  std::optional<std::string> transmission;
  std::optional<double> site_altitude;
  std::optional<std::uint64_t> seed;
  std::optional<std::string> profile;
  std::optional<WavelengthRange> wavelengths;
  std::vector<std::pair<std::string, std::string>> envelope_files;
  OutputMode mode = OutputMode::thinned;
  bool store_p = false;
  std::optional<std::string> scattering;
  std::optional<std::uint32_t> scattering_points;
  std::string input = "-";
  std::string output = "-";

  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      const auto value = [&](const char *name) -> std::string {
        if (i + 1 >= argc) {
          throw UsageError(std::string(name) + " needs a value");
        }
        return argv[++i];
      };
      if (arg == "--transmission") {
        transmission = value("--transmission");
      } else if (arg == "--site-altitude") {
        site_altitude = parse_double(value("--site-altitude"), "--site-altitude");
      } else if (arg == "--profile") {
        profile = value("--profile");
      } else if (arg == "--envelope") {
        const std::string spec = value("--envelope");
        const auto equals = spec.find('=');
        if (equals == std::string::npos || equals == 0 || equals + 1 == spec.size()) {
          throw UsageError("--envelope needs NAME=FILE: " + spec);
        }
        envelope_files.emplace_back(spec.substr(0, equals), spec.substr(equals + 1));
      } else if (arg == "--wavelength-range") {
        const double lo = parse_double(value("--wavelength-range"), "--wavelength-range LO");
        if (i + 1 >= argc) {
          throw UsageError("--wavelength-range needs LO and HI");
        }
        const double hi = parse_double(argv[++i], "--wavelength-range HI");
        wavelengths = WavelengthRange{lo, hi};
      } else if (arg == "--seed") {
        seed = parse_seed(value("--seed"));
      } else if (arg == "--mode") {
        const std::string text = value("--mode");
        if (text == "thinned") {
          mode = OutputMode::thinned;
        } else if (text == "weighted") {
          mode = OutputMode::weighted;
        } else {
          throw UsageError("--mode is thinned or weighted: " + text);
        }
      } else if (arg == "--store-p") {
        store_p = true;
      } else if (arg == "--scattering") {
        scattering = value("--scattering");
      } else if (arg == "--scattering-points") {
        scattering_points = parse_points(value("--scattering-points"));
      } else if (arg == "-i") {
        input = value("-i");
      } else if (arg == "-o") {
        output = value("-o");
      } else if (arg == "-h" || arg == "--help") {
        usage();
        return 0;
      } else {
        throw UsageError("unknown argument: " + arg);
      }
    }
    if (!transmission || !site_altitude || !seed) {
      throw UsageError("--transmission, --site-altitude and --seed are required");
    }
    if (scattering_points && !scattering) {
      throw UsageError("--scattering-points needs --scattering");
    }

    PropagatorConfig config;
    config.transmission.emplace(TransmissionTable::from_file(*transmission, *site_altitude));
    config.transmission_source = *transmission;
    config.site_altitude_m = *site_altitude;
    if (profile) {
      config.profile.emplace(AtmosphereProfile::from_file(*profile));
      config.profile_source = *profile;
    }
    for (const auto &[name, path] : envelope_files) {
      config.envelope.push_back(read_envelope_file(name, path));
    }
    config.wavelength_range = wavelengths;
    config.seed = *seed;
    config.mode = mode;
    config.store_p = store_p;
    if (scattering) {
      config.scattering.enabled = true;
      config.scattering.points_per_photon = scattering_points.value_or(8);
      config.scattering_table.emplace(ScatteringTable::from_file(*scattering));
      config.scattering_source = *scattering;
    }

    std::ifstream in_file;
    if (input != "-") {
      in_file.open(input, std::ios::binary);
      if (!in_file) {
        throw std::runtime_error("cannot open input stream: " + input);
      }
    }
    std::ofstream out_file;
    if (output != "-") {
      out_file.open(output, std::ios::binary);
      if (!out_file) {
        throw std::runtime_error("cannot open output stream: " + output);
      }
    }
    std::istream &in =
        input == "-" ? static_cast<std::istream &>(std::cin) : static_cast<std::istream &>(in_file);
    std::ostream &out = output == "-" ? static_cast<std::ostream &>(std::cout)
                                      : static_cast<std::ostream &>(out_file);
    const PropagationSummary summary = propagate_stream(in, out, std::move(config));
    out.flush();
    std::fprintf(stderr, "dlufd-propagate: photons in %llu, out %llu (weight %.6g), dropped %llu\n",
                 static_cast<unsigned long long>(summary.photons_in),
                 static_cast<unsigned long long>(summary.photons_out), summary.weight_out,
                 static_cast<unsigned long long>(summary.dropped));
  } catch (const UsageError &error) {
    std::fprintf(stderr, "dlufd-propagate: %s\n", error.what());
    usage();
    return 2;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "dlufd-propagate: error: %s\n", error.what());
    return 1;
  }
  return 0;
}
