// corsika2ps: CORSIKA IACT EventIO (file or stdin) to a photon stream on stdout.

#include "gstrickd/simtel.hpp"
#include "gstrickd/stdin_buffer.hpp"
#include "gstrickd/version.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr std::size_t buffer_size = std::size_t{1} << 20;

void usage() {
  std::fprintf(stderr, "usage: corsika2ps [input|-]  (writes the stream to stdout)\n");
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
  std::string input = "-";
  bool have_input = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    }
    if (arg.size() > 1 && arg[0] == '-') {
      std::fprintf(stderr, "corsika2ps: unknown option %s\n", arg.c_str());
      usage();
      return 2;
    }
    if (have_input) {
      std::fprintf(stderr, "corsika2ps: only one input is accepted\n");
      usage();
      return 2;
    }
    input = arg;
    have_input = true;
  }

  try {
    static std::vector<char> out_buffer(buffer_size);
    std::cout.rdbuf()->pubsetbuf(out_buffer.data(),
                                 static_cast<std::streamsize>(out_buffer.size()));
    gstrickd::StdinBuffer stdin_buffer;
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
    const std::string producer = "corsika2ps " + std::to_string(gstrickd::version_major) + "." +
                                 std::to_string(gstrickd::version_minor) + "." +
                                 std::to_string(gstrickd::version_patch);
    gstrickd::corsika_to_stream(*in, std::cout, producer);
    std::cout.flush();
    if (!std::cout) {
      throw std::runtime_error("write error on stdout");
    }
  } catch (const std::exception &error) {
    std::fprintf(stderr, "corsika2ps: %s\n", error.what());
    return 1;
  }
  return 0;
}
