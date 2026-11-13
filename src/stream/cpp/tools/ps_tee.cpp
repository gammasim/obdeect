// obdeect stream tee: copy a photon stream from stdin to stdout and to each named file,
// record by record, flushing after every record. Photon columns are not decoded.

#include "obdeect/stream/stdin_buffer.hpp"
#include "obdeect/stream/stream.hpp"
#include "tool_version.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr std::size_t buffer_size = std::size_t{1} << 20;

// Checks the record order of a raw copy: one stream record first, stream_end
// last, nothing after stream_end. Records are checked before they are written.
class StreamOrder {
public:
  void accept(obdeect::stream::RecordType type) {
    if (ended_) {
      throw std::runtime_error("data after stream_end");
    }
    if (type == obdeect::stream::RecordType::stream) {
      if (header_seen_) {
        throw std::runtime_error("second stream record");
      }
      header_seen_ = true;
      return;
    }
    if (!header_seen_) {
      throw std::runtime_error("stream does not start with a stream record");
    }
    ended_ = type == obdeect::stream::RecordType::stream_end;
  }

  // Called at EOF on a record boundary.
  void finish() const {
    if (!header_seen_) {
      throw std::runtime_error("empty input (no stream record)");
    }
    if (!ended_) {
      throw std::runtime_error("no stream_end (truncated stream)");
    }
  }

private:
  bool header_seen_{false};
  bool ended_{false};
};

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  std::vector<std::string> names;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--version") {
      print_version("obdeect stream tee");
      return 0;
    }
    if (arg == "-h" || arg == "--help") {
      std::fprintf(
          stderr, "usage: obdeect stream tee file...  (copies stdin to stdout and to each file)\n");
      return 0;
    }
    if (arg.size() > 1 && arg[0] == '-') {
      std::fprintf(stderr, "obdeect stream tee: unknown option %s\n", arg.c_str());
      return 2;
    }
    names.push_back(arg);
  }

  try {
    static std::vector<char> out_buffer(buffer_size);
    obdeect::stream::StdinBuffer stdin_buffer;
    std::istream in(&stdin_buffer);
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
    StreamOrder order;
    while (obdeect::stream::read_record_raw(in, record)) {
      const auto type = static_cast<obdeect::stream::RecordType>(
          obdeect::stream::stream_detail::load_le(record.data() + 4, 2));
      order.accept(type);
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
    order.finish();
  } catch (const std::exception &error) {
    std::fprintf(stderr, "obdeect stream tee: %s\n", error.what());
    return 1;
  }
  return 0;
}
