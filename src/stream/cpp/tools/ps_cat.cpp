// ps-cat: concatenate photon streams into one valid stream on stdout, or print
// a summary of the streams.

#include "gstrickd/stream.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr std::size_t buffer_size = std::size_t{1} << 20;

struct Summary {
  std::array<std::uint64_t, 16> records{};
  std::map<std::uint32_t, std::uint64_t> photons;
  std::map<std::uint32_t, double> weights;
};

void usage() { std::fprintf(stderr, "usage: ps-cat [--summary] [file|-]...\n"); }

// Copies the records of one stream (checked: framing, format major, stream_end).
void copy_stream(std::istream &in, const std::string &name, gstrickd::Writer &writer) {
  gstrickd::Reader reader(in);
  gstrickd::Record record;
  while (reader.next(record)) {
    switch (record.type) {
    case gstrickd::RecordType::stream:
      break; // one stream header per output
    case gstrickd::RecordType::photons:
      writer.write_chunk(record.chunk);
      break;
    case gstrickd::RecordType::blob:
      writer.write_blob(record.blob_eventio_type, record.blob_eventio_id, record.blob);
      break;
    default:
      writer.write_metadata(record.type, record.metadata);
      break;
    }
  }
  if (!reader.ended()) {
    throw std::runtime_error(name + ": no stream_end (truncated stream)");
  }
  if (in.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error(name + ": data after stream_end");
  }
}

void summarize_stream(std::istream &in, const std::string &name, Summary &summary) {
  gstrickd::Reader reader(in);
  gstrickd::Record record;
  while (reader.next(record)) {
    summary.records[static_cast<std::size_t>(record.type)] += 1;
    if (record.type == gstrickd::RecordType::photons) {
      const auto id = record.chunk.telescope_id;
      summary.photons[id] += record.chunk.rows();
      double sum = 0.0;
      for (const float w : record.chunk[gstrickd::Column::w]) {
        sum += w;
      }
      summary.weights[id] += sum;
    }
  }
  if (!reader.ended()) {
    throw std::runtime_error(name + ": no stream_end (truncated stream)");
  }
  if (in.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error(name + ": data after stream_end");
  }
}

// Calls fn(istream&, name) for each input; "-" is stdin.
template <class Fn> void for_each_input(const std::vector<std::string> &inputs, Fn &&fn) {
  for (const auto &name : inputs) {
    if (name == "-") {
      fn(std::cin, name);
      continue;
    }
    std::vector<char> buffer(buffer_size);
    std::ifstream file;
    file.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    file.open(name, std::ios::binary);
    if (!file) {
      throw std::runtime_error("cannot open " + name);
    }
    fn(file, name);
  }
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
  bool summary = false;
  std::vector<std::string> inputs;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--summary") {
      summary = true;
    } else if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    } else if (arg.size() > 1 && arg[0] == '-') {
      std::fprintf(stderr, "ps-cat: unknown option %s\n", arg.c_str());
      usage();
      return 2;
    } else {
      inputs.push_back(arg);
    }
  }
  if (inputs.empty()) {
    inputs.emplace_back("-");
  }

  try {
    if (summary) {
      Summary result;
      for_each_input(inputs, [&](std::istream &in, const std::string &name) {
        summarize_stream(in, name, result);
      });
      std::printf("records:\n");
      for (std::size_t t = 0; t < result.records.size(); ++t) {
        if (result.records[t] > 0) {
          std::printf("  %-14s %llu\n", gstrickd::record_name(static_cast<gstrickd::RecordType>(t)),
                      static_cast<unsigned long long>(result.records[t]));
        }
      }
      std::printf("photons per telescope:\n");
      std::uint64_t total_photons = 0;
      double total_weight = 0.0;
      for (const auto &[id, count] : result.photons) {
        std::printf("  telescope %u: %llu photons, weight sum %.10g\n", id,
                    static_cast<unsigned long long>(count), result.weights[id]);
        total_photons += count;
        total_weight += result.weights[id];
      }
      std::printf("total: %llu photons, weight sum %.10g\n",
                  static_cast<unsigned long long>(total_photons), total_weight);
      return 0;
    }

    static std::vector<char> out_buffer(buffer_size);
    std::cout.rdbuf()->pubsetbuf(out_buffer.data(),
                                 static_cast<std::streamsize>(out_buffer.size()));
    gstrickd::Writer writer(std::cout, "ps-cat");
    for_each_input(
        inputs, [&](std::istream &in, const std::string &name) { copy_stream(in, name, writer); });
    writer.finish();
    std::cout.flush();
    if (!std::cout) {
      throw std::runtime_error("write error on stdout");
    }
  } catch (const std::exception &error) {
    std::fprintf(stderr, "ps-cat: %s\n", error.what());
    return 1;
  }
  return 0;
}
