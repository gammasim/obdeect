// ps-tee: copy a photon stream from stdin to stdout and to each named file,
// record by record, flushing after every record. Photon columns are not decoded.

#include "gstrickd/stream.hpp"

#include <cstdint>
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
  std::vector<std::string> names;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      std::fprintf(stderr, "usage: ps-tee file...  (copies stdin to stdout and to each file)\n");
      return 0;
    }
    if (arg.size() > 1 && arg[0] == '-') {
      std::fprintf(stderr, "ps-tee: unknown option %s\n", arg.c_str());
      return 2;
    }
    names.push_back(arg);
  }

  try {
    static std::vector<char> in_buffer(buffer_size);
    static std::vector<char> out_buffer(buffer_size);
    std::cin.rdbuf()->pubsetbuf(in_buffer.data(), static_cast<std::streamsize>(in_buffer.size()));
    std::cout.rdbuf()->pubsetbuf(out_buffer.data(),
                                 static_cast<std::streamsize>(out_buffer.size()));

    std::vector<std::vector<char>> buffers(names.size(), std::vector<char>(buffer_size));
    std::vector<std::ofstream> outputs(names.size());
    for (std::size_t i = 0; i < names.size(); ++i) {
      outputs[i].rdbuf()->pubsetbuf(buffers[i].data(),
                                    static_cast<std::streamsize>(buffers[i].size()));
      outputs[i].open(names[i], std::ios::binary);
      if (!outputs[i]) {
        throw std::runtime_error("cannot open " + names[i] + " for writing");
      }
    }

    std::vector<std::uint8_t> record;
    while (gstrickd::read_record_raw(std::cin, record)) {
      const auto *data = reinterpret_cast<const char *>(record.data());
      const auto size = static_cast<std::streamsize>(record.size());
      std::cout.write(data, size);
      std::cout.flush();
      for (auto &output : outputs) {
        output.write(data, size);
        output.flush();
      }
      if (!std::cout) {
        throw std::runtime_error("write error on stdout");
      }
      for (std::size_t i = 0; i < outputs.size(); ++i) {
        if (!outputs[i]) {
          throw std::runtime_error("write error on " + names[i]);
        }
      }
    }
  } catch (const std::exception &error) {
    std::fprintf(stderr, "ps-tee: %s\n", error.what());
    return 1;
  }
  return 0;
}
