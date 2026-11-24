// obdeect stream to-grisu: ground-frame photon stream (file or stdin) to GrISU/GrOptics photon text
// on stdout or -o FILE.

#include "obdeect/stream/grisu.hpp"
#include "obdeect/stream/parse.hpp"
#include "tool_io.hpp"
#include "tool_version.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

void usage() {
  std::fprintf(stderr, "usage: obdeect stream to-grisu [--telescope-config FILE] [--qe FACTOR] "
                       "[--expand-bunches] [--event-details] [-o FILE|-] [input|-]\n");
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  return run_stream_tool("obdeect stream to-grisu", argc, argv, usage, [&] {
    obdeect::stream::cli::Args args(argc, argv);
    std::string input = "-";
    std::string output = "-";
    std::string config_path;
    obdeect::stream::GrisuOptions options;
    options.producer = std::string("obdeect stream to-grisu ") + OBDEECT_STREAM_PACKAGE_VERSION;
    options.warnings = &std::cerr;
    bool have_input = false;
    std::string arg;
    while (args.next(arg)) {
      if (arg == "-h" || arg == "--help") {
        usage();
        return 0;
      }
      if (arg == "--expand-bunches") {
        options.expand_bunches = true;
      } else if (arg == "--event-details") {
        options.event_details = true;
      } else if (arg == "-o") {
        output = args.value(arg);
      } else if (arg == "--telescope-config") {
        config_path = args.value(arg);
      } else if (arg == "--qe") {
        options.qe = obdeect::stream::cli::number(args.value(arg), arg);
      } else if (arg.size() > 1 && arg[0] == '-') {
        throw obdeect::stream::cli::UsageError("unknown option " + arg);
      } else if (have_input) {
        throw obdeect::stream::cli::UsageError("only one input is accepted");
      } else {
        input = arg;
        have_input = true;
      }
    }

    if (!config_path.empty()) {
      std::ifstream config(config_path);
      if (!config) {
        throw std::runtime_error("cannot open telescope config " + config_path);
      }
      options.telescopes = obdeect::stream::read_tlloc(config);
    }

    obdeect::stream::tool::Input in{input};
    obdeect::stream::tool::Output out{output};
    obdeect::stream::stream_to_grisu(in.stream(), out.stream(), std::move(options));
    out.finish();
    return 0;
  });
}
