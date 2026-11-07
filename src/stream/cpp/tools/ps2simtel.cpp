// ps2simtel: photon stream (file or stdin) to sim_telarray input (EventIO) on stdout or -o FILE.

#include "gstrickd/simtel.hpp"
#include "gstrickd/stdin_buffer.hpp"

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
  std::fprintf(stderr, "usage: ps2simtel [--camera-map FILE] [-o FILE|-] [input|-]\n");
}

void set_binary_stdio() {
#ifdef _WIN32
  _setmode(_fileno(stdin), _O_BINARY);
  _setmode(_fileno(stdout), _O_BINARY);
#endif
}

// Prints the sim_telarray Bypass_Optics value for each converted telescope.
void report_bypass(const gstrickd::CameraMap &camera) {
  for (const auto &[id, telescope] : camera) {
    std::fprintf(stderr,
                 "ps2simtel: telescope %u: camera frame, mirror class %d, Bypass_Optics=%d, "
                 "camera rotation %g deg\n",
                 id, telescope.mirror_class, gstrickd::bypass_optics_level(telescope.mirror_class),
                 telescope.rotation_deg);
  }
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  set_binary_stdio();
  std::string input = "-";
  std::string output = "-";
  std::string camera_path;
  bool have_input = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    }
    if (arg == "--camera-map" || arg == "-o") {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "ps2simtel: %s needs an argument\n", arg.c_str());
        usage();
        return 2;
      }
      (arg == "-o" ? output : camera_path) = argv[++i];
      continue;
    }
    if (arg.size() > 1 && arg[0] == '-') {
      std::fprintf(stderr, "ps2simtel: unknown option %s\n", arg.c_str());
      usage();
      return 2;
    }
    if (have_input) {
      std::fprintf(stderr, "ps2simtel: only one input is accepted\n");
      usage();
      return 2;
    }
    input = arg;
    have_input = true;
  }

  try {
    gstrickd::CameraMap camera;
    if (!camera_path.empty()) {
      std::ifstream map_file(camera_path);
      if (!map_file) {
        throw std::runtime_error("cannot open camera map " + camera_path);
      }
      camera = gstrickd::read_camera_map(map_file);
      report_bypass(camera);
    }

    static std::vector<char> in_buffer(buffer_size);
    static std::vector<char> out_buffer(buffer_size);
    gstrickd::StdinBuffer stdin_buffer;
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

    gstrickd::stream_to_simtel(*in, *out, std::move(camera));
    out->flush();
    if (!*out) {
      throw std::runtime_error("write error on output");
    }
    if (out_file.is_open()) {
      out_file.close();
    }
  } catch (const std::exception &error) {
    std::fprintf(stderr, "ps2simtel: %s\n", error.what());
    return 1;
  }
  return 0;
}
