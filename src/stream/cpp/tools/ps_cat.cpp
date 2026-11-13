// ps-cat: concatenate photon streams into one valid stream on stdout, or print
// a summary of the streams.

#include "obdeect/stream/stdin_buffer.hpp"
#include "obdeect/stream/stream.hpp"
#include "tool_version.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr std::size_t buffer_size = std::size_t{1} << 20;

struct Summary {
  std::array<std::uint64_t, 16> records{};
  std::map<std::uint32_t, std::uint64_t> photons;
  std::map<std::uint32_t, double> weights;
};

void usage() { std::fprintf(stderr, "usage: ps-cat [--summary] [file|-]...\n"); }

// Copies the records of one stream (checked: framing, format major, stream_end).
void copy_stream(std::istream &in, const std::string &name, obdeect::stream::Writer &writer) {
  obdeect::stream::Reader reader(in);
  obdeect::stream::Record record;
  while (reader.next(record)) {
    switch (record.type) {
    case obdeect::stream::RecordType::stream:
      break; // one stream header per output
    case obdeect::stream::RecordType::photons:
      writer.write_chunk(record.chunk);
      break;
    case obdeect::stream::RecordType::blob:
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
  if (!reader.at_clean_end()) {
    throw std::runtime_error(name + ": data after stream_end");
  }
}

void summarize_stream(std::istream &in, const std::string &name, Summary &summary) {
  obdeect::stream::Reader reader(in);
  obdeect::stream::Record record;
  while (reader.next(record)) {
    summary.records[static_cast<std::size_t>(record.type)] += 1;
    if (record.type == obdeect::stream::RecordType::photons) {
      const auto id = record.chunk.telescope_id;
      summary.photons[id] += record.chunk.rows();
      double sum = 0.0;
      for (const float w : record.chunk[obdeect::stream::Column::w]) {
        sum += w;
      }
      summary.weights[id] += sum;
    }
  }
  if (!reader.ended()) {
    throw std::runtime_error(name + ": no stream_end (truncated stream)");
  }
  if (!reader.at_clean_end()) {
    throw std::runtime_error(name + ": data after stream_end");
  }
}

// Calls fn(istream&, name) for each input; "-" is stdin.
template <class Fn>
void for_each_input(const std::vector<std::string> &inputs, std::istream &stdin_stream, Fn &&fn) {
  for (const auto &name : inputs) {
    if (name == "-") {
      fn(stdin_stream, name);
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

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  bool summary = false;
  std::vector<std::string> inputs;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--summary") {
      summary = true;
    } else if (arg == "--version") {
      print_version("ps-cat");
      return 0;
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
    obdeect::stream::StdinBuffer stdin_buffer;
    std::istream stdin_stream(&stdin_buffer);
    if (summary) {
      Summary result;
      for_each_input(inputs, stdin_stream, [&](std::istream &in, const std::string &name) {
        summarize_stream(in, name, result);
      });
      std::printf("records:\n");
      for (std::size_t t = 0; t < result.records.size(); ++t) {
        if (result.records[t] > 0) {
          std::printf("  %-14s %llu\n", obdeect::stream::record_name(static_cast<obdeect::stream::RecordType>(t)),
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
    obdeect::stream::Writer writer(std::cout, "ps-cat");
    for_each_input(inputs, stdin_stream, [&](std::istream &in, const std::string &name) {
      copy_stream(in, name, writer);
    });
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
