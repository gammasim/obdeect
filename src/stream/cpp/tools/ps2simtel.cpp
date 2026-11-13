// obdeect stream to-simtel: photon stream (file or stdin) to sim_telarray input (EventIO) on stdout
// or -o FILE.

#include "obdeect/stream/simtel.hpp"
#include "obdeect/stream/stdin_buffer.hpp"
#include "tool_version.hpp"

#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr std::size_t buffer_size = std::size_t{1} << 20;

void usage() {
  std::fprintf(stderr, "usage: obdeect stream to-simtel [--camera-map FILE] [--atmprof FILE] "
                       "[--allow-pending FACTOR]... [--scaled-factor FACTOR]... [-o FILE|-] "
                       "[input|-]\n");
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  return run_stream_tool("obdeect stream to-simtel", argc, argv, usage, [&] {
    obdeect::stream::cli::Args args(argc, argv);
    std::string input = "-";
    std::string output = "-";
    std::string camera_path;
    std::string atmprof_path;
    std::set<std::string> allow_pending;
    std::set<std::string> scaled;
    bool have_input = false;
    std::string arg;
    while (args.next(arg)) {
      if (arg == "-h" || arg == "--help") {
        usage();
        return 0;
      }
      if (arg == "-o") {
        output = args.value(arg);
      } else if (arg == "--allow-pending") {
        allow_pending.insert(args.value(arg));
      } else if (arg == "--scaled-factor") {
        scaled.insert(args.value(arg));
      } else if (arg == "--atmprof") {
        atmprof_path = args.value(arg);
      } else if (arg == "--camera-map") {
        camera_path = args.value(arg);
      } else if (arg.size() > 1 && arg[0] == '-') {
        throw obdeect::stream::cli::UsageError("unknown option " + arg);
      } else if (have_input) {
        throw obdeect::stream::cli::UsageError("only one input is accepted");
      } else {
        input = arg;
        have_input = true;
      }
    }

    obdeect::stream::SimtelOptions options;
    options.warnings = &std::cerr;
    options.allow_pending = std::move(allow_pending);
    options.scaled = std::move(scaled);
    if (!camera_path.empty()) {
      std::ifstream map_file(camera_path);
      if (!map_file) {
        throw std::runtime_error("cannot open camera map " + camera_path);
      }
      options.camera = obdeect::stream::read_camera_map(map_file);
    }
    if (!atmprof_path.empty()) {
      std::ifstream table_file(atmprof_path);
      if (!table_file) {
        throw std::runtime_error("cannot open atmospheric table " + atmprof_path);
      }
      options.atmprof = obdeect::stream::read_atmprof_text(table_file, atmprof_path);
    }

    static std::vector<char> in_buffer(buffer_size);
    static std::vector<char> out_buffer(buffer_size);
    obdeect::stream::StdinBuffer stdin_buffer;
    std::istream stdin_stream(&stdin_buffer);
    std::ifstream in_file;
    std::istream *in = &stdin_stream;
    if (input != "-") {
      in_file.rdbuf()->pubsetbuf(in_buffer.data(), static_cast<std::streamsize>(in_buffer.size()));
      in_file.open(input, std::ios::binary);
      if (!in_file) {
        throw std::runtime_error("cannot open " + input);
      }
      in = &in_file;
    }

    std::ofstream out_file;
    std::ostream *out = &std::cout;
    if (output == "-") {
      std::cout.rdbuf()->pubsetbuf(out_buffer.data(),
                                   static_cast<std::streamsize>(out_buffer.size()));
    } else {
      out_file.rdbuf()->pubsetbuf(out_buffer.data(),
                                  static_cast<std::streamsize>(out_buffer.size()));
      out_file.open(output, std::ios::binary);
      if (!out_file) {
        throw std::runtime_error("cannot open " + output + " for writing");
      }
      out = &out_file;
    }

    obdeect::stream::stream_to_simtel(*in, *out, std::move(options));
    out->flush();
    if (!*out) {
      throw std::runtime_error("write error on output");
    }
    if (out_file.is_open()) {
      out_file.close();
    }
    return 0;
  });
}
