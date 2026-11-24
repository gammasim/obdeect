// obdeect stream cat: concatenate photon streams into one valid stream on stdout, or print
// a summary of the streams.

#include "obdeect/stream/stream.hpp"
#include "tool_io.hpp"
#include "tool_version.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Summary {
  std::array<std::uint64_t, 16> records{};
  std::map<std::uint32_t, std::uint64_t> photons;
  std::map<std::uint32_t, double> weights;
};

void usage() { std::fprintf(stderr, "usage: obdeect stream cat [--summary] [file|-]...\n"); }

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

} // namespace

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  return run_stream_tool("obdeect stream cat", argc, argv, usage, [&] {
    obdeect::stream::cli::Args args(argc, argv);
    bool summary = false;
    std::vector<std::string> inputs;
    std::string arg;
    while (args.next(arg)) {
      if (arg == "--summary") {
        summary = true;
      } else if (arg == "-h" || arg == "--help") {
        usage();
        return 0;
      } else if (arg.size() > 1 && arg[0] == '-') {
        throw obdeect::stream::cli::UsageError("unknown option " + arg);
      } else {
        inputs.push_back(arg);
      }
    }

    if (inputs.empty()) {
      inputs.emplace_back("-");
    }
    if (summary) {
      Summary result;
      for (const std::string &name : inputs) {
        obdeect::stream::tool::Input in{name};
        summarize_stream(in.stream(), name, result);
      }
      std::printf("records:\n");
      for (std::size_t t = 0; t < result.records.size(); ++t) {
        if (result.records[t] > 0) {
          std::printf("  %-14s %llu\n",
                      obdeect::stream::record_name(static_cast<obdeect::stream::RecordType>(t)),
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
      if (std::fflush(stdout) != 0 || std::ferror(stdout) != 0) {
        throw std::runtime_error("write error on standard output");
      }
      return 0;
    }

    obdeect::stream::tool::Output out{"-"};
    obdeect::stream::Writer writer(out.stream(), "obdeect stream cat");
    for (const std::string &name : inputs) {
      obdeect::stream::tool::Input in{name};
      copy_stream(in.stream(), name, writer);
    }
    writer.finish();
    out.finish();
    return 0;
  });
}
