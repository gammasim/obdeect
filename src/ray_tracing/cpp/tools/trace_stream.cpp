// obdeect ray-tracing trace-stream: trace a photon stream (G or T photons) to camera-frame (C)
// photons (docs/ray-tracing/photon-stream.md). Reads a stream from a file or stdin, writes the
// traced stream to a file or stdout, and prints a summary to stderr. Same options and output as
// obdeect.ray_tracing.stream.trace_stream (Python).

#include "obdeect/ray_tracing/stream_trace.hpp"
#include "obdeect/stream/cli.hpp"
#include "obdeect/stream/parse.hpp"
#include "obdeect/stream/stdin_buffer.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

namespace rt = obdeect::ray_tracing;
using obdeect::stream::cli::UsageError;

void usage() {
  std::fprintf(stderr,
               "usage: obdeect ray-tracing trace-stream --model [ID=]PATH [--model ID=PATH]...\n"
               "                       --seed N [--threads N] [--telescopes ID[,ID...]]\n"
               "                       [--pointing-zenith-deg DEG --pointing-azimuth-deg DEG]\n"
               "                       [--telescope-transmission full|angular]\n"
               "                       [-i FILE|-] [-o FILE|-]\n");
}

constexpr std::uint32_t max_telescope_id = 65535;

std::uint32_t telescope_id(long value) {
  if (value < 1 || value > static_cast<long>(max_telescope_id))
    throw rt::ModelError("telescope id " + std::to_string(value) + " outside 1.." +
                         std::to_string(max_telescope_id));
  return static_cast<std::uint32_t>(value);
}

// --telescopes IDS: comma-separated telescope ids; empty items are ignored.
std::set<std::uint32_t> parse_telescopes(const std::string &text) {
  std::set<std::uint32_t> ids;
  std::size_t start = 0;
  while (start <= text.size()) {
    const std::size_t comma = std::min(text.find(',', start), text.size());
    std::string item = text.substr(start, comma - start);
    item.erase(0, item.find_first_not_of(" \t"));
    item.erase(item.find_last_not_of(" \t") + 1);
    if (!item.empty()) {
      long id = 0;
      if (!obdeect::stream::parse::integer(item, id))
        throw UsageError("telescope ids must be integers: " + text);
      ids.insert(telescope_id(id));
    }
    start = comma + 1;
  }
  if (ids.empty())
    throw UsageError("give at least one telescope id");
  return ids;
}

// --model PATH (every telescope) or ID=PATH (one telescope).
std::pair<std::optional<std::uint32_t>, std::string> parse_model(const std::string &text) {
  const auto equals = text.find('=');
  if (equals == std::string::npos)
    return {std::nullopt, text};
  long id = 0;
  if (!obdeect::stream::parse::integer(text.substr(0, equals), id))
    throw UsageError("--model takes PATH or TELESCOPE_ID=PATH: " + text);
  return {telescope_id(id), text.substr(equals + 1)};
}

// --threads N: 0 uses every core, as in the Python command.
unsigned parse_threads(const std::string &text) {
  long value = 0;
  if (!obdeect::stream::parse::integer(text, value) || value < 0 || value > 256)
    throw UsageError("--threads needs an integer from 0 to 256: " + text);
  if (value == 0)
    return std::max(1U, std::thread::hardware_concurrency());
  return static_cast<unsigned>(value);
}

std::string read_file(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot open optical model: " + path);
  return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

bool same_file(const std::string &input, const std::string &output) {
  std::error_code error;
  const bool same = std::filesystem::equivalent(input, output, error);
  return !error && same;
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  return obdeect::stream::cli::run_tool(
      "obdeect ray-tracing trace-stream", obdeect::stream::cli::core_version(), argc, argv, usage,
      [&] {
        std::vector<std::pair<std::optional<std::uint32_t>, std::string>> models;
        std::optional<std::uint64_t> seed;
        unsigned threads = std::max(1U, std::thread::hardware_concurrency());
        std::optional<std::set<std::uint32_t>> telescopes;
        std::optional<double> zenith, azimuth;
        auto part = rt::stream_trace_detail::TransmissionPart::full;
        std::string input = "-";
        std::string output = "-";
        obdeect::stream::cli::Args args(argc, argv);
        std::string arg;
        while (args.next(arg)) {
          const auto value = [&](const char *name) -> std::string { return args.value(name); };
          if (arg == "--model") {
            models.push_back(parse_model(value("--model")));
          } else if (arg == "--seed") {
            std::uint64_t parsed = 0;
            const std::string text = value("--seed");
            if (!obdeect::stream::parse::integer(text, parsed))
              throw UsageError("--seed needs a non-negative integer: " + text);
            seed = parsed;
          } else if (arg == "--threads") {
            threads = parse_threads(value("--threads"));
          } else if (arg == "--telescopes") {
            telescopes = parse_telescopes(value("--telescopes"));
          } else if (arg == "--pointing-zenith-deg") {
            zenith = obdeect::stream::cli::number(value("--pointing-zenith-deg"),
                                                  "--pointing-zenith-deg");
          } else if (arg == "--pointing-azimuth-deg") {
            azimuth = obdeect::stream::cli::number(value("--pointing-azimuth-deg"),
                                                   "--pointing-azimuth-deg");
          } else if (arg == "--telescope-transmission") {
            const std::string text = value("--telescope-transmission");
            if (text == "full")
              part = rt::stream_trace_detail::TransmissionPart::full;
            else if (text == "angular")
              part = rt::stream_trace_detail::TransmissionPart::angular;
            else
              throw UsageError("--telescope-transmission is full or angular: " + text);
          } else if (arg == "-i" || arg == "--input") {
            input = value("--input");
          } else if (arg == "-o" || arg == "--output") {
            output = value("--output");
          } else if (arg == "-h" || arg == "--help") {
            usage();
            return 0;
          } else {
            throw UsageError("unknown argument: " + arg);
          }
        }
        if (models.empty() || !seed)
          throw UsageError("--model and --seed are required");
        if (zenith.has_value() != azimuth.has_value())
          throw rt::ModelError(
              "give both --pointing-zenith-deg and --pointing-azimuth-deg, or neither");

        rt::StreamTraceSettings settings;
        settings.seed = *seed;
        settings.threads = threads;
        settings.part = part;
        settings.version = obdeect::stream::cli::core_version();
        if (zenith)
          settings.pointing = std::pair{*zenith, *azimuth};

        // One model for every telescope, or the model of each telescope id (files loaded once).
        std::shared_ptr<const rt::stream_trace_detail::Optics> shared;
        std::map<std::uint32_t, std::shared_ptr<const rt::stream_trace_detail::Optics>> per;
        if (models.size() == 1 && !models.front().first) {
          shared = rt::stream_trace_detail::make_optics(
              rt::stream_trace_detail::load_stream_model(read_file(models.front().second)), part,
              settings.version);
          settings.telescopes = telescopes;
        } else {
          std::map<std::filesystem::path, std::shared_ptr<const rt::stream_trace_detail::Optics>>
              loaded;
          for (const auto &[id, path] : models) {
            if (!id)
              throw rt::ModelError(
                  "give one --model PATH, or --model TELESCOPE_ID=PATH for each telescope");
            if (per.count(*id))
              throw rt::ModelError("a telescope id is given in more than one --model");
            const auto key = std::filesystem::weakly_canonical(path);
            if (!loaded.count(key))
              loaded[key] = rt::stream_trace_detail::make_optics(
                  rt::stream_trace_detail::load_stream_model(read_file(path)), part,
                  settings.version);
            per[*id] = loaded[key];
          }
          std::set<std::uint32_t> selected;
          if (telescopes) {
            selected = *telescopes;
          } else {
            for (const auto &[id, optics] : per)
              selected.insert(id);
          }
          std::string missing;
          for (const auto id : selected)
            if (!per.count(id))
              missing += (missing.empty() ? "" : ", ") + std::to_string(id);
          if (!missing.empty())
            throw rt::ModelError("telescopes without an optical model: [" + missing + "]");
          for (auto it = per.begin(); it != per.end();)
            it = selected.count(it->first) ? std::next(it) : per.erase(it);
          settings.telescopes = selected;
        }

        if (input != "-" && output != "-" && same_file(input, output))
          throw rt::ModelError("input and output are the same file: " + output);
        std::ifstream in_file;
        if (input != "-") {
          in_file.open(input, std::ios::binary);
          if (!in_file)
            throw std::runtime_error("cannot open input stream: " + input);
        }
        std::ofstream out_file;
        if (output != "-") {
          out_file.open(output, std::ios::binary);
          if (!out_file)
            throw std::runtime_error("cannot open output stream: " + output);
        }
        // std::cin is read one character at a time by libc++ once sync_with_stdio is off.
        obdeect::stream::StdinBuffer stdin_buffer;
        std::istream stdin_stream(&stdin_buffer);
        std::istream &in = input == "-" ? stdin_stream : static_cast<std::istream &>(in_file);
        std::ostream &out = output == "-" ? static_cast<std::ostream &>(std::cout)
                                          : static_cast<std::ostream &>(out_file);
        const auto counts = rt::trace_stream(in, out, std::move(shared), std::move(per), settings);
        out.flush();
        if (!out)
          throw std::runtime_error("write error on the output stream");
        // Standard output may carry the stream, so the summary goes to standard error.
        std::fprintf(stderr, "trace-stream: %llu events, %llu photons in, %llu out\n",
                     static_cast<unsigned long long>(counts.events),
                     static_cast<unsigned long long>(counts.photons_in),
                     static_cast<unsigned long long>(counts.photons_out));
        return 0;
      });
}
