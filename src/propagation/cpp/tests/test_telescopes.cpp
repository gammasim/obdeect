#include "check.hpp"

#include "obdeect/propagation/propagator.hpp"
#include "obdeect/propagation/scattering.hpp"
#include "obdeect/propagation/transmission.hpp"

#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef OBDEECT_PROPAGATION_TEST_DATA_DIR
#error "OBDEECT_PROPAGATION_TEST_DATA_DIR must be defined by CMake"
#endif

// Telescope selection (propagation run --telescopes, PropagatorConfig::telescopes): the output
// for the selected telescopes equals the output of the full run with the other blocks removed,
// and photons of other telescopes are not decoded.

namespace {

using obdeect::stream::Column;
using obdeect::stream::Frame;
using obdeect::stream::RecordType;

constexpr double site_m = 2000.0;
constexpr std::size_t telescope_count = 3;

std::string make_stream(const std::function<void(obdeect::stream::Writer &)> &body) {
  std::ostringstream out;
  {
    obdeect::stream::Writer writer(out, "test");
    body(writer);
    writer.finish();
  }
  return out.str();
}

// Stream of two events, each with telescopes 1..3 (G photons, n per telescope). Telescope 2 has a
// telescope-level efficiency status (qe applied upstream) after its telescope record.
std::string input_stream(std::size_t n) {
  return make_stream([&](obdeect::stream::Writer &writer) {
    obdeect::stream::Metadata run;
    run["run_number"] = std::int64_t{1};
    run["observation_level_m"] = site_m;
    run["telescope_x_m"] = std::vector<double>{100.0, -50.0, 20.0};
    run["telescope_y_m"] = std::vector<double>{0.0, 300.0, -120.0};
    run["telescope_z_m"] = std::vector<double>{0.0, 0.0, 0.0};
    run["telescope_r_m"] = std::vector<double>{12.0, 12.0, 12.0};
    run["wavelength_min_nm"] = 250.0;
    run["wavelength_max_nm"] = 700.0;
    writer.write_metadata(RecordType::run, run);
    for (std::int64_t event = 7; event < 9; ++event) {
      obdeect::stream::Metadata event_keys;
      event_keys["event_number"] = event;
      writer.write_metadata(RecordType::event, event_keys);
      for (std::uint32_t telescope = 1; telescope <= telescope_count; ++telescope) {
        obdeect::stream::Metadata open;
        open["telescope_id"] = static_cast<std::int64_t>(telescope);
        open["frame"] = static_cast<std::int64_t>(Frame::ground);
        writer.write_metadata(RecordType::telescope, open);
        if (telescope == 2) {
          obdeect::stream::Metadata status;
          status["factors"] = std::string("qe");
          status["qe.state"] = std::string("applied");
          status["qe.source"] = std::string("upstream");
          writer.write_metadata(RecordType::efficiency_status, status);
        }
        for (std::size_t first = 0; first < n; first += 50) {
          const std::size_t rows = std::min<std::size_t>(50, n - first);
          obdeect::stream::PhotonChunk chunk;
          chunk.telescope_id = telescope;
          chunk.frame = Frame::ground;
          chunk.first_index = first;
          chunk.resize(rows, obdeect::stream::column_bit(Column::x) |
                                 obdeect::stream::column_bit(Column::y) |
                                 obdeect::stream::column_bit(Column::cx) |
                                 obdeect::stream::column_bit(Column::cy) |
                                 obdeect::stream::column_bit(Column::t) |
                                 obdeect::stream::column_bit(Column::w) |
                                 obdeect::stream::column_bit(Column::zem));
          for (std::size_t r = 0; r < rows; ++r) {
            const auto i = static_cast<double>(first + r + 1000 * telescope + 100000 * event);
            chunk[Column::x][r] = static_cast<float>(300.0 * std::sin(0.37 * i));
            chunk[Column::y][r] = static_cast<float>(300.0 * std::cos(0.23 * i));
            chunk[Column::cx][r] = static_cast<float>(0.4 * std::sin(0.001 * i));
            chunk[Column::cy][r] = static_cast<float>(0.4 * std::cos(0.0013 * i));
            chunk[Column::t][r] = static_cast<float>(std::fmod(i, 13.0));
            chunk[Column::zem][r] = static_cast<float>(2000.0 + 1000.0 * std::fmod(i, 7.0));
            chunk[Column::w][r] = static_cast<float>(1.0 + std::fmod(i, 4.0));
          }
          writer.write_chunk(chunk);
        }
        obdeect::stream::Metadata end;
        end["telescope_id"] = static_cast<std::int64_t>(telescope);
        end["photons"] = 0.0;
        writer.write_metadata(RecordType::telescope_end, end);
      }
      writer.write_metadata(RecordType::event_end, {});
    }
    writer.write_metadata(RecordType::run_end, {});
  });
}

std::vector<obdeect::stream::Record> records_of(const std::string &bytes) {
  std::istringstream in(bytes);
  obdeect::stream::Reader reader(in);
  std::vector<obdeect::stream::Record> records;
  obdeect::stream::Record record;
  while (reader.next(record)) {
    records.push_back(record);
  }
  CHECK(reader.ended());
  return records;
}

// The records written again with the writer of propagate_stream (same producer): equal records
// give equal bytes.
std::string reserialize(const std::vector<obdeect::stream::Record> &records) {
  std::ostringstream out;
  {
    obdeect::stream::Writer writer(out, "obdeect propagation");
    for (const auto &record : records) {
      switch (record.type) {
      case RecordType::stream:
        break;
      case RecordType::photons:
        writer.write_chunk(record.chunk);
        break;
      case RecordType::blob:
        writer.write_blob(record.blob_eventio_type, record.blob_eventio_id, record.blob);
        break;
      default:
        writer.write_metadata(record.type, record.metadata);
        break;
      }
    }
    writer.finish();
  }
  return out.str();
}

// The records of the full run that a selection keeps: run-level records, and the blocks (telescope,
// photons, efficiency status, telescope_end) of the selected telescopes.
std::vector<obdeect::stream::Record>
expected_records(const std::vector<obdeect::stream::Record> &full,
                 const std::vector<std::uint32_t> &selection) {
  const auto kept = [&](std::uint32_t id) {
    return std::find(selection.begin(), selection.end(), id) != selection.end();
  };
  std::vector<obdeect::stream::Record> out;
  bool in_block = false;
  std::uint32_t block = 0;
  for (const auto &record : full) {
    switch (record.type) {
    case RecordType::telescope:
      block =
          static_cast<std::uint32_t>(std::get<std::int64_t>(record.metadata.at("telescope_id")));
      in_block = true;
      if (kept(block)) {
        out.push_back(record);
      }
      break;
    case RecordType::telescope_end:
      if (kept(block)) {
        out.push_back(record);
      }
      in_block = false;
      break;
    case RecordType::efficiency_status:
      if (!in_block || kept(block)) {
        out.push_back(record);
      }
      break;
    case RecordType::photons:
      if (kept(record.chunk.telescope_id)) {
        out.push_back(record);
      }
      break;
    default:
      out.push_back(record);
      break;
    }
  }
  return out;
}

obdeect::propagation::PropagatorConfig config_for(const obdeect::propagation::TransmissionTable &t,
                                                  std::vector<std::uint32_t> selection = {}) {
  obdeect::propagation::PropagatorConfig config;
  config.transmission.emplace(t);
  config.transmission_source = "test table";
  config.site_altitude_m = site_m;
  config.seed = 20261107;
  if (!selection.empty()) {
    config.telescopes = std::move(selection);
  }
  return config;
}

std::string propagate(const std::string &input, obdeect::propagation::PropagatorConfig config,
                      obdeect::propagation::PropagationSummary *summary = nullptr) {
  std::istringstream in(input);
  std::ostringstream out;
  const auto result = obdeect::propagation::propagate_stream(in, out, std::move(config));
  if (summary != nullptr) {
    *summary = result;
  }
  return out.str();
}

template <class Function> bool throws_invalid(Function function) {
  try {
    function();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

obdeect::propagation::TransmissionTable table() {
  return obdeect::propagation::TransmissionTable::from_file(
      std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/atm_trans_2000_1_10_0_0_2000.dat", site_m);
}

} // namespace

int main() {
  return run_checks([] {
    const std::string input = input_stream(120);
    const auto full_bytes = propagate(input, config_for(table()));
    const auto full = records_of(full_bytes);

    // One telescope: the output is the full output with the other blocks removed, byte for byte.
    {
      obdeect::propagation::PropagationSummary summary;
      const auto selected = propagate(input, config_for(table(), {2}), &summary);
      CHECK(selected == reserialize(expected_records(full, {2})));
      CHECK(selected.size() < full_bytes.size());
      // Counts are for the selected telescope only (input rows of telescope 2: 2 events x 120).
      CHECK(summary.photons_in == 240);
      CHECK(summary.photons_out > 0);
    }

    // Order and repeats in the selection do not matter: {3, 1, 3} keeps telescopes 1 and 3.
    {
      const auto selected = propagate(input, config_for(table(), {3, 1, 3}));
      CHECK(selected == reserialize(expected_records(full, {1, 3})));
    }

    // Keeping every telescope gives the unselected output.
    {
      const auto selected = propagate(input, config_for(table(), {1, 2, 3}));
      CHECK(selected == full_bytes);
    }

    // The stream reader with a selection returns photons of other telescopes without columns.
    {
      std::istringstream in(input);
      obdeect::stream::Reader reader(in);
      reader.keep_photons_of({2});
      obdeect::stream::Record record;
      std::size_t skipped = 0;
      std::size_t decoded = 0;
      while (reader.next(record)) {
        if (record.type != RecordType::photons) {
          continue;
        }
        if (record.chunk.telescope_id == 2) {
          decoded += record.chunk.rows();
        } else {
          CHECK(record.chunk.rows() == 0);
          skipped += 1;
        }
      }
      CHECK(decoded == 240);
      CHECK(skipped > 0);
    }

    // Invalid selections and combinations are rejected.
    CHECK(throws_invalid([&] { propagate(input, config_for(table(), {0})); }));
    CHECK(throws_invalid([&] { propagate(input, config_for(table(), {4})); }));
    {
      auto config = config_for(table());
      config.telescopes = std::vector<std::uint32_t>{};
      CHECK(throws_invalid([&] { propagate(input, config); }));
    }
    {
      auto config = config_for(table(), {2});
      config.scattering.enabled = true;
      config.scattering_table.emplace(obdeect::propagation::ScatteringTable::from_file(
          std::string(OBDEECT_PROPAGATION_TEST_DATA_DIR) + "/scattering_exponential.dat"));
      CHECK(throws_invalid([&] { propagate(input, config); }));
    }
  });
}
