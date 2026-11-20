#pragma once

// Frame G to frame sky: the exact inverse of sky_to_ground (sky.hpp) for photons of a CORSIKA
// stream, used by `obdeect propagation to-sky` (docs/propagation/ground-to-sky.md).
//
// Per row, with the telescope position (X_t, Y_t, Z_t) [m], the observation level h_obs [m a.s.l.],
// h_t = h_obs + Z_t, cz_m = sqrt(1 - cx^2 - cy^2) (photons move downwards, cz = -cz_m) and
// zenith = atan2(sqrt(cx^2 + cy^2), cz_m):
//   s     = (zem - h_t) / cz_m                         straight path from the emission to the plane
//   x_sky = X_t + x - s cx,  y_sky = Y_t + y - s cy,  z_sky = zem - h_obs
//   t_sky = t - (s / c + delay),
//           delay = straight_path_air_delay_ns(profile, h_t, zem, zenith) with a profile, else 0
// The stage (sky_to_ground with refraction_applied = 1) adds exactly these terms to the arrival
// time, so `to-sky` followed by `run --no-extinction` (same profile option) reproduces x, y, cx,
// cy, t and zem to float32 rounding. Bending of the photons is not inverted: a run with
// refraction_applied = 0 and a profile is rejected. The emission point is the straight
// back-projection along the arrival direction.

#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/timing.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::propagation {

// Geometry of one telescope for frame G photons. Heights are m a.s.l. (observation_m) or m above
// the observation level (telescope_z_m).
struct GroundPlane {
  double telescope_x_m{};
  double telescope_y_m{};
  double telescope_z_m{};
  double observation_m{};
  const AtmosphereProfile *profile{}; // straight-path air delay; null: vacuum
};

// Converts one row. Input: G position x, y, direction cx, cy, emission altitude zem [m a.s.l.] and
// time t [ns]. Output: sky x, y, z (above the observation level) and t. Returns false if the row
// cannot be converted: a non-finite input, zem <= h_t, or a horizontal photon (cx^2 + cy^2 = 1).
// Throws std::invalid_argument if cx^2 + cy^2 > 1. No allocation.
inline bool ground_to_sky_row(const GroundPlane &plane, double x, double y, double cx, double cy,
                              double zem, double t, double &x_out, double &y_out, double &z_out,
                              double &t_out) {
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(cx) || !std::isfinite(cy) ||
      !std::isfinite(zem) || !std::isfinite(t)) {
    return false;
  }
  if (!obdeect::stream::direction_in_unit_disc(cx, cy)) {
    throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
  }
  const double r2 = cx * cx + cy * cy;
  const double cz_magnitude = std::sqrt(std::max(0.0, 1.0 - r2));
  const double h_t = plane.observation_m + plane.telescope_z_m;
  if (!(cz_magnitude > 0.0) || !(zem > h_t)) {
    return false;
  }
  const double path = (zem - h_t) / cz_magnitude;
  x_out = plane.telescope_x_m + x - path * cx;
  y_out = plane.telescope_y_m + y - path * cy;
  z_out = zem - plane.observation_m;
  double delay_ns = path / obdeect::stream::speed_of_light_m_per_ns;
  if (plane.profile != nullptr) {
    const double zenith = std::atan2(std::sqrt(r2), cz_magnitude);
    delay_ns += straight_path_air_delay_ns(*plane.profile, h_t, zem, zenith);
  }
  t_out = t - delay_ns;
  return true;
}

// Rows [begin, end) of a chunk as a chunk of their own: first_index continues the row index.
inline obdeect::stream::PhotonChunk ground_slice(const obdeect::stream::PhotonChunk &chunk,
                                                 std::size_t begin, std::size_t end) {
  obdeect::stream::PhotonChunk piece;
  piece.telescope_id = chunk.telescope_id;
  piece.frame = chunk.frame;
  piece.first_index = chunk.first_index + begin;
  const std::uint32_t mask = chunk.mask();
  piece.resize(end - begin, mask);
  for (std::uint32_t i = 0; i < obdeect::stream::n_columns; ++i) {
    if ((mask >> i) & 1U) {
      const auto &source = chunk.columns[i];
      std::copy(source.begin() + static_cast<std::ptrdiff_t>(begin),
                source.begin() + static_cast<std::ptrdiff_t>(end), piece.columns[i].begin());
    }
  }
  return piece;
}

// Converts the rows of a frame G chunk to frame sky. Rows that cannot be converted (see
// ground_to_sky_row, or cz present and not downwards, or non-finite cz) are removed and added to
// dropped. The result is one sky chunk per run of consecutive converted rows, each with the
// first_index of its first row, so that the row indices of the photons (and thus the random draws
// of the later stages) do not change. Columns: x, y, t converted; z added; zem, w, wavelength, p
// and cz are kept when present. Input: x y cx cy t zem w (cz optional). Throws
// std::invalid_argument for frames other than G, a missing zem column, an invalid plane or an
// inconsistent cz column.
inline std::vector<obdeect::stream::PhotonChunk>
ground_to_sky(const obdeect::stream::PhotonChunk &input, const GroundPlane &plane,
              std::uint64_t &dropped) {
  using obdeect::stream::Column;
  if (input.frame != obdeect::stream::Frame::ground) {
    throw std::invalid_argument("ground_to_sky needs photons in frame G");
  }
  input.validate();
  if (!std::isfinite(plane.telescope_x_m) || !std::isfinite(plane.telescope_y_m) ||
      !std::isfinite(plane.telescope_z_m) || !std::isfinite(plane.observation_m)) {
    throw std::invalid_argument("telescope position and observation level must be finite");
  }
  const std::size_t rows = input.rows();
  if (rows == 0) {
    return {};
  }
  if (!input.has(Column::zem)) {
    throw std::invalid_argument("G-frame photons need the zem column (emission altitude)");
  }
  const bool has_cz = input.has(Column::cz);
  const auto &x_in = input[Column::x];
  const auto &y_in = input[Column::y];
  const auto &cx_in = input[Column::cx];
  const auto &cy_in = input[Column::cy];
  const auto &t_in = input[Column::t];
  const auto &zem_in = input[Column::zem];
  const auto &cz_in = input[Column::cz];

  obdeect::stream::PhotonChunk converted = input;
  converted.frame = obdeect::stream::Frame::sky;
  converted[Column::z].assign(rows, 0.0F);
  auto &x_out = converted[Column::x];
  auto &y_out = converted[Column::y];
  auto &z_out = converted[Column::z];
  auto &t_out = converted[Column::t];
  std::vector<std::uint8_t> keep(rows, 0);

  for (std::size_t row = 0; row < rows; ++row) {
    bool ok = true;
    if (has_cz) {
      const double cz = cz_in[row];
      if (!std::isfinite(cz) || cz >= 0.0) {
        ok = false; // moving upwards or horizontally, or not finite
      } else if (std::isfinite(cx_in[row]) && std::isfinite(cy_in[row]) &&
                 obdeect::stream::direction_in_unit_disc(cx_in[row], cy_in[row]) &&
                 !obdeect::stream::cz_matches_direction(cz, cx_in[row], cy_in[row])) {
        throw std::invalid_argument("cz column does not match cx, cy (cz = -sqrt(1-cx^2-cy^2))");
      }
    }
    if (!ok) {
      continue;
    }
    double x_sky = 0.0;
    double y_sky = 0.0;
    double z_sky = 0.0;
    double t_sky = 0.0;
    if (ground_to_sky_row(plane, x_in[row], y_in[row], cx_in[row], cy_in[row], zem_in[row],
                          t_in[row], x_sky, y_sky, z_sky, t_sky)) {
      keep[row] = 1;
      x_out[row] = static_cast<float>(x_sky);
      y_out[row] = static_cast<float>(y_sky);
      z_out[row] = static_cast<float>(z_sky);
      t_out[row] = static_cast<float>(t_sky);
    }
  }

  std::vector<obdeect::stream::PhotonChunk> out;
  std::uint64_t kept = 0;
  std::size_t row = 0;
  while (row < rows) {
    if (keep[row] == 0) {
      ++row;
      continue;
    }
    std::size_t end = row;
    while (end < rows && keep[end] != 0) {
      ++end;
    }
    out.push_back(ground_slice(converted, row, end));
    kept += end - row;
    row = end;
  }
  dropped += static_cast<std::uint64_t>(rows) - kept;
  return out;
}

struct ToSkySummary {
  std::uint64_t photons_in{};  // rows of the input photon records
  std::uint64_t photons_out{}; // rows written (sky frame)
  std::uint64_t dropped{};     // rows that cannot be converted (see ground_to_sky)
};

namespace to_sky_detail {

inline const obdeect::stream::Value &field(const obdeect::stream::Metadata &metadata,
                                           const char *key) {
  const auto it = metadata.find(key);
  if (it == metadata.end()) {
    throw std::invalid_argument(std::string("stream lacks the key '") + key + "'");
  }
  return it->second;
}

inline double number(const obdeect::stream::Metadata &metadata, const char *key) {
  const auto &value = field(metadata, key);
  if (const auto *real = std::get_if<double>(&value)) {
    return *real;
  }
  if (const auto *integer = std::get_if<std::int64_t>(&value)) {
    return static_cast<double>(*integer);
  }
  throw std::invalid_argument(std::string("key '") + key + "' is not a number");
}

inline std::int64_t integer(const obdeect::stream::Metadata &metadata, const char *key) {
  const auto &value = field(metadata, key);
  if (const auto *integer_value = std::get_if<std::int64_t>(&value)) {
    return *integer_value;
  }
  throw std::invalid_argument(std::string("key '") + key + "' is not an integer");
}

inline const std::vector<double> &array(const obdeect::stream::Metadata &metadata,
                                        const char *key) {
  const auto &value = field(metadata, key);
  if (const auto *values = std::get_if<std::vector<double>>(&value)) {
    return *values;
  }
  throw std::invalid_argument(std::string("key '") + key + "' is not a float array");
}

// Run geometry needed for the conversion (checked in run_geometry).
struct RunSite {
  double observation_m{};
  std::vector<double> x_m, y_m, z_m; // index = telescope_id - 1
};

inline RunSite run_geometry(const obdeect::stream::Metadata &metadata,
                            const std::optional<AtmosphereProfile> &profile) {
  if (metadata.count("curvature_applied") != 0) {
    const std::int64_t curvature = integer(metadata, "curvature_applied");
    if (curvature == 1) {
      throw std::invalid_argument("input photons have CORSIKA CURVED geometry "
                                  "(curvature_applied = 1): no flat sky geometry to convert to");
    }
    if (curvature != 0) {
      throw std::invalid_argument("curvature_applied must be 0 or 1");
    }
  }
  if (metadata.count("refraction_applied") != 0) {
    const std::int64_t refraction = integer(metadata, "refraction_applied");
    if (refraction != 0 && refraction != 1) {
      throw std::invalid_argument("refraction_applied must be 0 or 1");
    }
    if (profile && refraction != 1) {
      throw std::invalid_argument("--profile needs refraction_applied = 1 in the run (photons "
                                  "already refracted); refraction_applied = 0 is not inverted");
    }
  } else if (profile) {
    throw std::invalid_argument("run lacks refraction_applied (needed with a profile)");
  }
  RunSite site;
  site.observation_m = number(metadata, "observation_level_m");
  if (!std::isfinite(site.observation_m)) {
    throw std::invalid_argument("observation_level_m must be finite");
  }
  site.x_m = array(metadata, "telescope_x_m");
  site.y_m = array(metadata, "telescope_y_m");
  site.z_m = array(metadata, "telescope_z_m");
  if (site.y_m.size() != site.x_m.size() || site.z_m.size() != site.x_m.size()) {
    throw std::invalid_argument("telescope position arrays differ in length");
  }
  for (std::size_t i = 0; i < site.x_m.size(); ++i) {
    if (!std::isfinite(site.x_m[i]) || !std::isfinite(site.y_m[i]) || !std::isfinite(site.z_m[i])) {
      throw std::invalid_argument("telescope positions must be finite");
    }
  }
  return site;
}

} // namespace to_sky_detail

// Stream level: reads a photon stream, converts the frame G photons of every telescope to frame
// sky and writes the stream. Run, event, blob and other records pass through unchanged (run:
// telescope positions checked; telescope records: frame = sky). The telescope position is taken
// from the run keys telescope_x_m, telescope_y_m and telescope_z_m (the keys the atmosphere stage
// uses); telescope record keys x_m, y_m, z_m, if present, must agree with them. Throws
// std::invalid_argument for invalid input and std::runtime_error for a malformed or truncated
// stream.
inline ToSkySummary ground_to_sky_stream(std::istream &in, std::ostream &out,
                                         const std::optional<AtmosphereProfile> &profile,
                                         const std::string &producer = "obdeect propagation") {
  using obdeect::stream::Frame;
  using obdeect::stream::RecordType;
  using to_sky_detail::integer;
  constexpr double position_tolerance_m = 1e-6;

  obdeect::stream::Writer writer(out, producer);
  obdeect::stream::Reader reader(in);
  ToSkySummary summary;
  std::optional<to_sky_detail::RunSite> run;
  std::optional<std::uint32_t> open_telescope; // telescope whose photons may follow
  obdeect::stream::Record record;
  while (reader.next(record)) {
    switch (record.type) {
    case RecordType::stream:
    case RecordType::stream_end:
      break; // the writer writes its own stream and stream_end records
    case RecordType::run:
      run = to_sky_detail::run_geometry(record.metadata, profile);
      writer.write_metadata(record.type, record.metadata);
      break;
    case RecordType::telescope: {
      if (!run) {
        throw std::invalid_argument("telescope record before run record");
      }
      const std::int64_t id = integer(record.metadata, "telescope_id");
      if (id < 1 || static_cast<std::size_t>(id) > run->x_m.size()) {
        throw std::invalid_argument("telescope_id outside the run's telescope positions (1-based)");
      }
      const std::int64_t frame = integer(record.metadata, "frame");
      if (frame != static_cast<std::int64_t>(Frame::ground)) {
        throw std::invalid_argument("to-sky needs photons in frame G (telescope frame " +
                                    std::to_string(frame) + ")");
      }
      const auto index = static_cast<std::size_t>(id - 1);
      const char *keys[] = {"x_m", "y_m", "z_m"};
      const double expected[] = {run->x_m[index], run->y_m[index], run->z_m[index]};
      for (std::size_t k = 0; k < 3; ++k) {
        if (record.metadata.count(keys[k]) != 0 &&
            std::abs(to_sky_detail::number(record.metadata, keys[k]) - expected[k]) >
                position_tolerance_m) {
          throw std::invalid_argument(std::string("telescope record ") + keys[k] +
                                      " differs from the run's telescope position");
        }
      }
      open_telescope = static_cast<std::uint32_t>(id);
      obdeect::stream::Metadata out_metadata = record.metadata;
      out_metadata["frame"] = static_cast<std::int64_t>(Frame::sky);
      writer.write_metadata(record.type, out_metadata);
      break;
    }
    case RecordType::photons: {
      if (!run || !open_telescope || record.chunk.telescope_id != *open_telescope) {
        throw std::invalid_argument("photons of a telescope that is not open");
      }
      const auto &chunk = record.chunk;
      summary.photons_in += chunk.rows();
      if (chunk.rows() == 0) {
        break;
      }
      const auto index = static_cast<std::size_t>(*open_telescope - 1);
      const GroundPlane plane{run->x_m[index], run->y_m[index], run->z_m[index], run->observation_m,
                              profile ? &*profile : nullptr};
      for (const auto &piece : ground_to_sky(chunk, plane, summary.dropped)) {
        writer.write_chunk(piece);
        summary.photons_out += piece.rows();
      }
      break;
    }
    case RecordType::telescope_end: {
      open_telescope.reset();
      writer.write_metadata(record.type, record.metadata);
      break;
    }
    case RecordType::blob:
      writer.write_blob(record.blob_eventio_type, record.blob_eventio_id, record.blob);
      break;
    default:
      writer.write_metadata(record.type, record.metadata);
      break;
    }
  }
  if (!reader.ended()) {
    throw std::runtime_error("truncated photon stream (no stream_end record)");
  }
  writer.finish();
  return summary;
}

} // namespace obdeect::propagation
