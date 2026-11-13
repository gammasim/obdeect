// obdeect stream from-corsika: CORSIKA IACT EventIO (file or stdin) to a photon stream on stdout.

#include "obdeect/stream/simtel.hpp"
#include "obdeect/stream/stdin_buffer.hpp"
#include "obdeect/stream/version.hpp"
#include "tool_version.hpp"

#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr std::size_t buffer_size = std::size_t{1} << 20;

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
      if (arg.size() > 1 && arg[0] == '-')
        throw obdeect::stream::cli::UsageError("unknown option " + arg);
      if (have_input)
        throw obdeect::stream::cli::UsageError("only one input is accepted");
      input = arg;
      have_input = true;
    }

    static std::vector<char> out_buffer(buffer_size);
    std::cout.rdbuf()->pubsetbuf(out_buffer.data(),
                                 static_cast<std::streamsize>(out_buffer.size()));
    obdeect::stream::StdinBuffer stdin_buffer;
    std::istream stdin_stream(&stdin_buffer);
    static std::vector<char> in_buffer(buffer_size);
    std::ifstream file;
    std::istream *in = &stdin_stream;
    if (input != "-") {
      file.rdbuf()->pubsetbuf(in_buffer.data(), static_cast<std::streamsize>(in_buffer.size()));
      file.open(input, std::ios::binary);
      if (!file) {
        throw std::runtime_error("cannot open " + input);
      }
      in = &file;
    }
    const std::string producer = "obdeect stream from-corsika " +
                                 std::to_string(obdeect::stream::version_major) + "." +
                                 std::to_string(obdeect::stream::version_minor) + "." +
                                 std::to_string(obdeect::stream::version_patch);
    obdeect::stream::corsika_to_stream(*in, std::cout, producer);
    std::cout.flush();
    if (!std::cout) {
      throw std::runtime_error("write error on stdout");
    }
    return 0;
  });
}
