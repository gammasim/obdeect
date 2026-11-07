#pragma once

// Early thinning with the detector envelope. A photon survives with probability
// p = p_atm(lambda) * E(lambda), where E is the envelope (the upper bound over the
// later wavelength-dependent factors, e.g. QE x max_qe_rel, see docs/extinction.md).
// sim_telarray applies the same cut before tracing (sim_telarray.c:2914).

#include "dlufd/streams.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/random.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace dlufd {

// Largest bunch weight handled per row. Above it the sampling would reuse Philox counters
// of the wavelength stream; split such bunches upstream.
inline constexpr std::uint64_t max_bunch_photons = 2048;

// Envelope E(lambda): piecewise linear in wavelength, values in [0, 1].
class Envelope {
public:
  Envelope(std::vector<double> wavelength_nm, std::vector<double> value)
      : wavelength_(std::move(wavelength_nm)), value_(std::move(value)) {
    if (wavelength_.empty() || wavelength_.size() != value_.size()) {
      throw std::invalid_argument("envelope needs equal, non-empty wavelength and value arrays");
    }
    for (std::size_t i = 0; i < wavelength_.size(); ++i) {
      if (!std::isfinite(wavelength_[i]) || !std::isfinite(value_[i]) || value_[i] < 0.0 ||
          value_[i] > 1.0) {
        throw std::invalid_argument("envelope values must be finite and within [0, 1]");
      }
      if (i > 0 && wavelength_[i] <= wavelength_[i - 1]) {
        throw std::invalid_argument("envelope wavelengths must be strictly increasing");
      }
    }
  }

  // Linear interpolation; outside the table the edge value is returned (callers check
  // coverage before use).
  [[nodiscard]] double operator()(double lambda_nm) const {
    if (lambda_nm <= wavelength_.front()) {
      return value_.front();
    }
    if (lambda_nm >= wavelength_.back()) {
      return value_.back();
    }
    const auto hi = static_cast<std::size_t>(
        std::upper_bound(wavelength_.begin(), wavelength_.end(), lambda_nm) - wavelength_.begin());
    const std::size_t lo = hi - 1;
    const double t = (lambda_nm - wavelength_[lo]) / (wavelength_[hi] - wavelength_[lo]);
    return value_[lo] * (1.0 - t) + value_[hi] * t;
  }

  [[nodiscard]] double min_wavelength_nm() const noexcept { return wavelength_.front(); }
  [[nodiscard]] double max_wavelength_nm() const noexcept { return wavelength_.back(); }

private:
  std::vector<double> wavelength_;
  std::vector<double> value_;
};

// Deviate number `index` of a photon key: uniform4 block index / 4, lane index % 4.
[[nodiscard]] inline double deviate(const gstrickd::DrawKey &key, std::uint64_t index) {
  return gstrickd::uniform4(key, static_cast<std::uint32_t>(index / 4))[index % 4];
}

// Number of photons leaving a bunch of weight w with acceptance p:
// Binomial(floor(w), p) + Bernoulli((w - floor(w)) p). The binomial is the sum of
// floor(w) independent Bernoulli(p) trials (exact for any p), trial i using deviate i.
// The fractional part uses deviate floor(w). Requires 0 <= w <= max_bunch_photons, 0 <= p <= 1.
[[nodiscard]] inline std::uint64_t surviving_photons(double w, double p,
                                                     const gstrickd::DrawKey &key) {
  const double whole = std::floor(w);
  const auto n = static_cast<std::uint64_t>(whole);
  std::uint64_t survivors = 0;
  for (std::uint64_t block = 0; block * 4 < n; ++block) {
    const auto u = gstrickd::uniform4(key, static_cast<std::uint32_t>(block));
    for (std::uint64_t lane = 0; lane < 4 && block * 4 + lane < n; ++lane) {
      survivors += u[lane] < p ? 1U : 0U;
    }
  }
  const double fraction = w - whole;
  if (fraction > 0.0 && deviate(key, n) < fraction * p) {
    ++survivors;
  }
  return survivors;
}

// Thinned mode. Every surviving photon becomes one row with w = 1 and the other columns
// copied from its bunch. p (acceptance, per input row) is written as column p if write_p.
// Output photons are numbered from out_first_index; the result is split into chunks of at
// most 2^16 rows. The output does not depend on how the input is split into chunks, as
// long as first_index of each input chunk is its offset in the telescope's photons.
[[nodiscard]] inline std::vector<gstrickd::PhotonChunk>
thin_chunk(const gstrickd::PhotonChunk &in, std::span<const double> p, std::uint64_t seed,
           std::uint32_t event, std::uint32_t telescope, std::uint64_t out_first_index,
           bool write_p) {
  using gstrickd::Column;
  if (in.has(Column::p)) {
    throw std::invalid_argument("input already has an acceptance column p");
  }
  if (p.size() != in.rows()) {
    throw std::invalid_argument("acceptance needs one value per row");
  }
  // Source row of every output photon.
  std::vector<std::size_t> source;
  for (std::size_t row = 0; row < in.rows(); ++row) {
    if (!(p[row] >= 0.0 && p[row] <= 1.0)) {
      throw std::invalid_argument("acceptance must be within [0, 1]");
    }
    const double w = in[Column::w][row];
    if (!(w >= 0.0 && w <= static_cast<double>(max_bunch_photons))) {
      throw std::invalid_argument("bunch weight must be within [0, max_bunch_photons]");
    }
    const gstrickd::DrawKey key{seed,      stage_atmosphere,     event,
                                telescope, in.first_index + row, stream_thinning};
    const std::uint64_t count = surviving_photons(w, p[row], key);
    source.insert(source.end(), count, row);
  }

  std::uint32_t mask = in.mask() | gstrickd::column_bit(Column::w);
  if (write_p) {
    mask |= gstrickd::column_bit(Column::p);
  }
  std::vector<gstrickd::PhotonChunk> out;
  for (std::size_t first = 0; first < source.size(); first += gstrickd::max_chunk_rows) {
    const std::size_t rows = std::min<std::size_t>(gstrickd::max_chunk_rows, source.size() - first);
    gstrickd::PhotonChunk chunk;
    chunk.telescope_id = in.telescope_id;
    chunk.frame = in.frame;
    chunk.first_index = out_first_index + first;
    chunk.resize(rows, mask);
    for (std::size_t j = 0; j < rows; ++j) {
      const std::size_t row = source[first + j];
      for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
        const auto column = static_cast<Column>(c);
        if (column == Column::w || column == Column::p || !in.has(column)) {
          continue;
        }
        chunk[column][j] = in[column][row];
      }
      chunk[Column::w][j] = 1.0F;
      if (write_p) {
        chunk[Column::p][j] = static_cast<float>(p[row]);
      }
    }
    out.push_back(std::move(chunk));
  }
  return out;
}

// Weighted mode. Each row is kept with w multiplied by p; p is written as column p if
// write_p. Returns the chunk in place; all input is checked before the chunk changes.
inline void weight_chunk(gstrickd::PhotonChunk &chunk, std::span<const double> p, bool write_p) {
  using gstrickd::Column;
  if (chunk.has(Column::p)) {
    throw std::invalid_argument("input already has an acceptance column p");
  }
  if (p.size() != chunk.rows()) {
    throw std::invalid_argument("acceptance needs one value per row");
  }
  auto &w = chunk[Column::w];
  for (std::size_t row = 0; row < chunk.rows(); ++row) {
    if (!(p[row] >= 0.0 && p[row] <= 1.0)) {
      throw std::invalid_argument("acceptance must be within [0, 1]");
    }
    if (!(w[row] >= 0.0F && std::isfinite(w[row]))) {
      throw std::invalid_argument("bunch weight must be finite and non-negative");
    }
  }
  for (std::size_t row = 0; row < chunk.rows(); ++row) {
    w[row] = static_cast<float>(static_cast<double>(w[row]) * p[row]);
  }
  if (write_p) {
    auto &column = chunk[Column::p];
    column.resize(chunk.rows());
    for (std::size_t row = 0; row < chunk.rows(); ++row) {
      column[row] = static_cast<float>(p[row]);
    }
  }
}

} // namespace dlufd
