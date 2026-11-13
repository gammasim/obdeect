// obdeect propagation run: the atmosphere stage on a photon stream (docs/propagation.md).
// Reads a stream from a file or stdin, writes the processed stream to a file or stdout, and
// prints a summary to stderr.

#include "obdeect/propagation/propagator.hpp"
#include "obdeect/stream/cli.hpp"
#include "obdeect/stream/parse.hpp"

using namespace obdeect::propagation;

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void usage() {
  std::fprintf(
      stderr,
      "usage: obdeect propagation run (--transmission FILE | --no-extinction) --site-altitude M\n"
      "                       [--profile FILE]\n"
      "                       [--envelope NAME=FILE]... [--wavelength-range LO HI]\n"
      "                       --seed N [--mode thinned|weighted] [--store-p]\n"
      "                       [--scattering FILE] [--scattering-points N]\n"
      "                       [--pathlen M] [--telescope-pathlen ID=M]...\n"
      "                       [-i FILE|-] [-o FILE|-]\n");
}

using obdeect::stream::cli::UsageError;

double parse_double(const std::string &text, const char *what) {
  return obdeect::stream::cli::number(text, what);
}

// Path inside a telescope [m]: finite and not negative.
double parse_pathlen(const std::string &text, const char *what) {
  const double value = parse_double(text, what);
  if (!std::isfinite(value) || value < 0.0) {
    throw UsageError(std::string(what) + " must be finite and not negative: " + text);
  }
  return value;
}

// ID=M for --telescope-pathlen: telescope_id 1-based.
std::pair<std::uint32_t, double> parse_telescope_pathlen(const std::string &spec) {
  const auto equals = spec.find('=');
  if (equals == std::string::npos || equals == 0 || equals + 1 == spec.size()) {
    throw UsageError("--telescope-pathlen needs ID=M: " + spec);
  }
  const std::string id_text = spec.substr(0, equals);
  long id = 0;
  if (!obdeect::stream::parse::integer(id_text, id) || id < 1 ||
      id > static_cast<long>(obdeect::stream::max_telescope)) {
    throw UsageError("--telescope-pathlen: telescope ID is an integer from 1: " + id_text);
  }
  return {static_cast<std::uint32_t>(id),
          parse_pathlen(spec.substr(equals + 1), "--telescope-pathlen M")};
}

std::uint32_t parse_points(const std::string &text) {
  long value = 0;
  if (!obdeect::stream::parse::integer(text, value) || value < 1 ||
      value > static_cast<long>(max_points_per_photon)) {
    throw UsageError("--scattering-points needs an integer from 1 to 1024: " + text);
  }
  return static_cast<std::uint32_t>(value);
}

std::uint64_t parse_seed(const std::string &text) {
  std::uint64_t value = 0;
  if (!obdeect::stream::parse::integer(text, value)) {
    throw UsageError("invalid seed: " + text);
  }
  return value;
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  return obdeect::stream::cli::run_tool(
      "obdeect propagation run", obdeect::stream::cli::core_version(), argc, argv, usage, [&] {
        std::optional<std::string> transmission;
        std::optional<double> site_altitude;
        std::optional<std::uint64_t> seed;
        std::optional<std::string> profile;
        std::optional<WavelengthRange> wavelengths;
        std::vector<std::pair<std::string, std::string>> envelope_files;
        std::optional<OutputMode> mode;
        bool no_extinction = false;
        bool store_p = false;
        std::optional<std::string> scattering;
        std::optional<std::uint32_t> scattering_points;
        double pathlen = 0.0;
        bool pathlen_given = false;
        std::map<std::uint32_t, double> telescope_pathlen;
        std::string input = "-";
        std::string output = "-";

        {
          obdeect::stream::cli::Args args(argc, argv);
          std::string arg;
          while (args.next(arg)) {
            const auto value = [&](const char *name) -> std::string { return args.value(name); };
            if (arg == "--no-extinction") {
              no_extinction = true;
            } else if (arg == "--transmission") {
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
              const double hi = parse_double(value("--wavelength-range"), "--wavelength-range HI");
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
            } else if (arg == "--pathlen") {
              pathlen = parse_pathlen(value("--pathlen"), "--pathlen");
              pathlen_given = true;
            } else if (arg == "--telescope-pathlen") {
              const auto [id, metres] = parse_telescope_pathlen(value("--telescope-pathlen"));
              if (telescope_pathlen.count(id) != 0) {
                throw UsageError("--telescope-pathlen given twice for telescope " +
                                 std::to_string(id));
              }
              telescope_pathlen[id] = metres;
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
          if (!site_altitude || !seed || (!transmission && !no_extinction)) {
            throw UsageError("--transmission (or --no-extinction), --site-altitude and --seed are "
                             "required");
          }
          if (no_extinction) {
            // Geometry only: nothing is thinned or weighted, so these options make no sense.
            if (transmission || pathlen_given || !telescope_pathlen.empty()) {
              throw UsageError(
                  "--no-extinction cannot be combined with --transmission, --pathlen or "
                  "--telescope-pathlen");
            }
            if (!envelope_files.empty()) {
              throw UsageError("--no-extinction cannot be combined with --envelope");
            }
            if (mode == OutputMode::thinned) {
              throw UsageError("--no-extinction cannot be combined with --mode thinned");
            }
            if (scattering || scattering_points) {
              throw UsageError("--no-extinction cannot be combined with --scattering or "
                               "--scattering-points");
            }
            if (store_p) {
              throw UsageError("--no-extinction cannot be combined with --store-p");
            }
          }
          if (scattering_points && !scattering) {
            throw UsageError("--scattering-points needs --scattering");
          }

          PropagatorConfig config;
          config.geometry_only = no_extinction;
          if (transmission) {
            config.transmission.emplace(
                TransmissionTable::from_file(*transmission, *site_altitude));
            config.transmission_source = *transmission;
          }
          config.site_altitude_m = *site_altitude;
          config.default_pathlen_m = pathlen;
          config.telescope_pathlen_m = telescope_pathlen;
          if (profile) {
            config.profile.emplace(AtmosphereProfile::from_file(*profile));
            config.profile_source = *profile;
          }
          for (const auto &[name, path] : envelope_files) {
            config.envelope.push_back(read_envelope_file(name, path));
          }
          config.wavelength_range = wavelengths;
          config.seed = *seed;
          config.mode = mode.value_or(no_extinction ? OutputMode::weighted : OutputMode::thinned);
          config.store_p = store_p;
          if (scattering) {
            config.scattering.enabled = true;
            config.scattering.points_per_photon = scattering_points.value_or(8);
            config.scattering_table.emplace(ScatteringTable::from_file(*scattering));
            config.scattering_source = *scattering;
          }

          if (input != "-" && output != "-" && obdeect::propagation::same_file(input, output)) {
            throw UsageError("input and output are the same file: " + input);
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
          std::istream &in = input == "-" ? static_cast<std::istream &>(std::cin)
                                          : static_cast<std::istream &>(in_file);
          std::ostream &out = output == "-" ? static_cast<std::ostream &>(std::cout)
                                            : static_cast<std::ostream &>(out_file);
          const PropagationSummary summary = propagate_stream(in, out, std::move(config));
          out.flush();
          std::fprintf(
              stderr,
              "obdeect propagation run: photons in %llu, out %llu (weight %.6g), dropped %llu\n",
              static_cast<unsigned long long>(summary.photons_in),
              static_cast<unsigned long long>(summary.photons_out), summary.weight_out,
              static_cast<unsigned long long>(summary.dropped));
          return 0;
        }
      });
}
