// obdeect stream tee: copy a photon stream from stdin to stdout and to each named file,
// record by record, flushing after every record. Photon columns are not decoded.

#include "obdeect/stream/stream.hpp"
#include "tool_io.hpp"
#include "tool_version.hpp"

#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

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

void usage() {
  std::fprintf(stderr,
               "usage: obdeect stream tee file...  (copies stdin to stdout and to each file)\n");
}

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  return run_stream_tool("obdeect stream tee", argc, argv, usage, [&] {
    obdeect::stream::cli::Args args(argc, argv);
    std::vector<std::string> names;
    std::string arg;
    while (args.next(arg)) {
      if (arg == "-h" || arg == "--help") {
        usage();
        return 0;
      }
      if (arg.size() > 1 && arg[0] == '-') {
        throw obdeect::stream::cli::UsageError("unknown option " + arg);
      }
      names.push_back(arg);
    }

    obdeect::stream::tool::Input in{"-"};
    obdeect::stream::tool::Output out{"-"};
    std::vector<std::unique_ptr<obdeect::stream::tool::Output>> files;
    for (const std::string &name : names) {
      files.push_back(std::make_unique<obdeect::stream::tool::Output>(name));
    }

    std::vector<std::uint8_t> record;
    StreamOrder order;
    while (obdeect::stream::read_record_raw(in.stream(), record)) {
      const auto type = static_cast<obdeect::stream::RecordType>(
          obdeect::stream::stream_detail::load_le(record.data() + 4, 2));
      order.accept(type);
      const auto *data = reinterpret_cast<const char *>(record.data());
      const auto size = static_cast<std::streamsize>(record.size());
      out.stream().write(data, size);
      out.stream().flush();
      if (!out.stream()) {
        throw std::runtime_error("write error on standard output");
      }
      for (std::size_t i = 0; i < files.size(); ++i) {
        files[i]->stream().write(data, size);
        files[i]->stream().flush();
        if (!files[i]->stream()) {
          throw std::runtime_error("write error on " + names[i]);
        }
      }
    }
    order.finish();
    out.finish();
    for (const auto &file : files) {
      file->finish();
    }
    return 0;
  });
}
