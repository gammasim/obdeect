#pragma once

// Early thinning with the detector envelope. A photon survives with probability
// p = p_atm(lambda) * E(lambda), where E is the envelope (the upper bound over the
// later wavelength-dependent factors, e.g. QE x max_qe_rel, see docs/propagation/extinction.md).

#include "obdeect/propagation/parallel.hpp"
#include "obdeect/propagation/streams.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace obdeect::propagation {

// Largest bunch weight handled per row. The thinning draws of a bunch cost one Philox block per
// four photons, so a larger bunch would make one row too expensive; split such bunches upstream.
// (The draws of different streams and rows never share a counter: DrawKey.stream separates them.)
inline constexpr std::uint64_t max_bunch_photons = 2048;

// Envelope E(lambda): a piecewise-linear factor or a product of such factors, in [0, 1].
class Envelope {
public:
  Envelope(std::vector<double> wavelength_nm, std::vector<double> value,
           std::vector<Envelope> factors = {})
      : wavelength_(std::move(wavelength_nm)), value_(std::move(value)),
        factors_(std::move(factors)) {
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

  // Interpolate each factor before multiplying; outside a factor the edge value is returned.
  // Callers check wavelength coverage before use.
  [[nodiscard]] double operator()(double lambda_nm) const {
    if (!factors_.empty()) {
      double product = 1.0;
      for (const auto &factor : factors_) {
        product *= factor(lambda_nm);
      }
      return product;
    }
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
  std::vector<Envelope> factors_;
};

// Deviate number `index` of a photon key: uniform4 block index / 4, lane index % 4.
[[nodiscard]] inline double deviate(const obdeect::stream::DrawKey &key, std::uint64_t index) {
  return obdeect::stream::uniform4(key, static_cast<std::uint32_t>(index / 4))[index % 4];
}

// Number of photons leaving a bunch of weight w with acceptance p:
// Binomial(floor(w), p) + Bernoulli((w - floor(w)) p). The binomial is the sum of
// floor(w) independent Bernoulli(p) trials (exact for any p), trial i using deviate i.
// The fractional part uses deviate floor(w). Requires 0 <= w <= max_bunch_photons, 0 <= p <= 1.
[[nodiscard]] inline std::uint64_t surviving_photons(double w, double p,
                                                     const obdeect::stream::DrawKey &key) {
  const double whole = std::floor(w);
  const auto n = static_cast<std::uint64_t>(whole);
  std::uint64_t survivors = 0;
  for (std::uint64_t block = 0; block * 4 < n; ++block) {
    const auto u = obdeect::stream::uniform4(key, static_cast<std::uint32_t>(block));
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
// long as first_index of each input chunk is its offset in the telescope's photons (the
// caller checks that). The survivor counts and the copy of the output rows run on the threads
// of `pool`.
[[nodiscard]] inline std::vector<obdeect::stream::PhotonChunk>
thin_chunk(const obdeect::stream::PhotonChunk &in, std::span<const double> p, std::uint64_t seed,
           std::uint32_t event, std::uint32_t telescope, std::uint64_t out_first_index,
           bool write_p, std::uint32_t stream, RowPool &pool) {
  using obdeect::stream::Column;
  if (in.has(Column::p)) {
    throw std::invalid_argument("input already has an acceptance column p");
  }
  if (p.size() != in.rows()) {
    throw std::invalid_argument("acceptance needs one value per row");
  }
  // Survivors of every row (rows are independent), then the source row of every output photon.
  const std::size_t rows = in.rows();
  std::vector<std::uint64_t> counts(rows, 0);
  const auto &weight = in[Column::w];
  pool.run(rows, [&](std::size_t begin, std::size_t end) {
    for (std::size_t row = begin; row < end; ++row) {
      if (!(p[row] >= 0.0 && p[row] <= 1.0)) {
        throw std::invalid_argument("acceptance must be within [0, 1]");
      }
      const double w = weight[row];
      if (!(w >= 0.0 && w <= static_cast<double>(max_bunch_photons))) {
        throw std::invalid_argument("bunch weight must be within [0, max_bunch_photons]");
      }
      const obdeect::stream::DrawKey key{seed,      stage_atmosphere,     event,
                                         telescope, in.first_index + row, stream};
      counts[row] = surviving_photons(w, p[row], key);
    }
  });
  std::size_t total = 0;
  for (const std::uint64_t count : counts) {
    total += count;
  }
  std::vector<std::size_t> source;
  source.reserve(total);
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::uint64_t k = 0; k < counts[row]; ++k) {
      source.push_back(row);
    }
  }

  std::uint32_t mask = in.mask() | obdeect::stream::column_bit(Column::w);
  if (write_p) {
    mask |= obdeect::stream::column_bit(Column::p);
  }
  std::vector<obdeect::stream::PhotonChunk> out;
  for (std::size_t first = 0; first < source.size(); first += obdeect::stream::max_chunk_rows) {
    const std::size_t rows =
        std::min<std::size_t>(obdeect::stream::max_chunk_rows, source.size() - first);
    obdeect::stream::PhotonChunk chunk;
    chunk.telescope_id = in.telescope_id;
    chunk.frame = in.frame;
    chunk.first_index = out_first_index + first;
    chunk.resize(rows, mask);
    pool.run(rows, [&](std::size_t begin, std::size_t end) {
      for (std::size_t j = begin; j < end; ++j) {
        const std::size_t row = source[first + j];
        for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c) {
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
    });
    out.push_back(std::move(chunk));
  }
  return out;
}

// thin_chunk with a pool for this call only (threads: at least 1).
[[nodiscard]] inline std::vector<obdeect::stream::PhotonChunk>
thin_chunk(const obdeect::stream::PhotonChunk &in, std::span<const double> p, std::uint64_t seed,
           std::uint32_t event, std::uint32_t telescope, std::uint64_t out_first_index,
           bool write_p, std::uint32_t stream = stream_thinning, unsigned threads = 1) {
  RowPool pool(threads);
  return thin_chunk(in, p, seed, event, telescope, out_first_index, write_p, stream, pool);
}

// Weighted mode. Each row is kept with w multiplied by p; p is written as column p if
// write_p. Returns the chunk in place; all input is checked before the chunk changes.
inline void weight_chunk(obdeect::stream::PhotonChunk &chunk, std::span<const double> p,
                         bool write_p) {
  using obdeect::stream::Column;
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

} // namespace obdeect::propagation
