// obdeect stream from-corsika: CORSIKA IACT EventIO (file or stdin) to a photon stream on stdout.

#include "obdeect/stream/cli.hpp"
#include "obdeect/stream/simtel.hpp"
#include "tool_io.hpp"
#include "tool_version.hpp"

#include <cstdio>
#include <string>

namespace {

void usage() {
  std::fprintf(stderr,
               "usage: obdeect stream from-corsika [input|-]  (writes the stream to stdout)\n");
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  return run_stream_tool("obdeect stream from-corsika", argc, argv, usage, [&] {
    obdeect::stream::cli::Args args(argc, argv);
    std::string input = "-";
    bool have_input = false;
    std::string arg;
    while (args.next(arg)) {
      if (arg == "-h" || arg == "--help") {
        usage();
        return 0;
      }
      if (arg.size() > 1 && arg[0] == '-') {
        throw obdeect::stream::cli::UsageError("unknown option " + arg);
      }
      if (have_input) {
        throw obdeect::stream::cli::UsageError("only one input is accepted");
      }
      input = arg;
      have_input = true;
    }

    obdeect::stream::tool::Input in{input};
    obdeect::stream::tool::Output out{"-"};
    obdeect::stream::corsika_to_stream(in.stream(), out.stream(),
                                       std::string("obdeect stream from-corsika ") +
                                           OBDEECT_STREAM_PACKAGE_VERSION);
    out.finish();
    return 0;
  });
}
