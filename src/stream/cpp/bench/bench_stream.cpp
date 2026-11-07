// Throughput of the photon-stream writer and reader through a file.
// Usage: bench_stream [photons] [file]   (defaults: 1e8 photons, temp file)

#include "gstrickd/stream.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

constexpr std::uint64_t chunk_rows = std::uint64_t{1} << 16;

double seconds_since(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

} // namespace

int main(int argc, char **argv) {
  const std::uint64_t photons = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 100000000ULL;
  const std::filesystem::path path =
      argc > 2 ? std::filesystem::path(argv[2])
               : std::filesystem::temp_directory_path() / "gstrickd_bench.ps";
  const std::uint64_t n_chunks = (photons + chunk_rows - 1) / chunk_rows;

  gstrickd::PhotonChunk chunk;
  chunk.telescope_id = 1;
  chunk.frame = gstrickd::Frame::telescope;
  chunk.resize(chunk_rows, gstrickd::required_columns |
                               gstrickd::column_bit(gstrickd::Column::zem) |
                               gstrickd::column_bit(gstrickd::Column::wavelength));
  for (auto &column : chunk.columns) {
    for (std::size_t r = 0; r < column.size(); ++r) {
      column[r] = static_cast<float>(r) * 0.001F;
    }
  }

  std::vector<char> buffer(std::size_t{1} << 20);
  double write_seconds = 0.0;
  double bytes = 0.0; // file size in bytes
  {
    std::ofstream out;
    out.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    out.open(path, std::ios::binary);
    if (!out) {
      std::fprintf(stderr, "cannot open %s\n", path.string().c_str());
      return EXIT_FAILURE;
    }
    const auto start = std::chrono::steady_clock::now();
    gstrickd::Writer writer(out, "bench_stream");
    for (std::uint64_t i = 0; i < n_chunks; ++i) {
      chunk.first_index = i * chunk_rows;
      writer.write_chunk(chunk);
    }
    writer.finish();
    out.close();
    write_seconds = seconds_since(start);
    bytes = static_cast<double>(std::filesystem::file_size(path));
  }

  std::uint64_t read_photons = 0;
  double weight_sum = 0.0;
  double read_seconds = 0.0;
  {
    std::ifstream in;
    in.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    in.open(path, std::ios::binary);
    const auto start = std::chrono::steady_clock::now();
    gstrickd::Reader reader(in);
    gstrickd::Record record;
    while (reader.next(record)) {
      if (record.type == gstrickd::RecordType::photons) {
        read_photons += record.chunk.rows();
        for (const float w : record.chunk[gstrickd::Column::w]) {
          weight_sum += w;
        }
      }
    }
    read_seconds = seconds_since(start);
  }
  std::filesystem::remove(path);

  const double mb = bytes / 1e6;
  std::printf("photons written: %llu (%llu chunks of %llu rows, 8 columns)\n",
              static_cast<unsigned long long>(n_chunks * chunk_rows),
              static_cast<unsigned long long>(n_chunks),
              static_cast<unsigned long long>(chunk_rows));
  std::printf("file size: %.1f MB\n", mb);
  std::printf("write: %.2f s, %.0f MB/s\n", write_seconds, mb / write_seconds);
  std::printf("read:  %.2f s, %.0f MB/s (read %llu photons, weight sum %.6g)\n", read_seconds,
              mb / read_seconds, static_cast<unsigned long long>(read_photons), weight_sum);
  return read_photons == n_chunks * chunk_rows ? EXIT_SUCCESS : EXIT_FAILURE;
}
