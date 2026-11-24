// obdeect stream to-simtel: photon stream (file or stdin) to sim_telarray input (EventIO) on stdout
// or -o FILE.

#include "obdeect/stream/simtel.hpp"
#include "tool_io.hpp"
#include "tool_version.hpp"

#include <cstdio>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>

namespace {

void usage() {
  std::fprintf(stderr, "usage: obdeect stream to-simtel [--camera-map FILE] [--atmprof FILE] "
                       "[--allow-pending FACTOR]... [--allow-pending-atmosphere] "
                       "[--scaled-factor FACTOR]... [--keep-telpos-height] "
                       "[-o FILE|-] [input|-]\n");
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
    bool keep_telpos_height = false;
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
      } else if (arg == "--allow-pending-atmosphere") {
        allow_pending.insert("atmosphere");
      } else if (arg == "--scaled-factor") {
        scaled.insert(args.value(arg));
      } else if (arg == "--keep-telpos-height") {
        keep_telpos_height = true;
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
    options.keep_telpos_height = keep_telpos_height;
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

    obdeect::stream::tool::Input in{input};
    obdeect::stream::tool::Output out{output};
    obdeect::stream::stream_to_simtel(in.stream(), out.stream(), std::move(options));
    out.finish();
    return 0;
  });
}
