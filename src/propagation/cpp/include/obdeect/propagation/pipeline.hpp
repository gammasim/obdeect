#pragma once

// One chunk of the atmosphere step for telescope photons in frame ground:
// wavelength assignment, atmospheric acceptance, envelope, then weighted or thinned output.

#include "obdeect/propagation/extinction.hpp"
#include "obdeect/propagation/parallel.hpp"
#include "obdeect/propagation/streams.hpp"
#include "obdeect/propagation/thinning.hpp"
#include "obdeect/propagation/transmission.hpp"
#include "obdeect/propagation/wavelength.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace obdeect::propagation {

struct ExtinctionSettings {
  const TransmissionTable *table{};
  const TelescopeDepths *depths{};
  const Envelope *envelope{};
  WavelengthRange wavelengths{};
  std::uint64_t seed{};
  std::uint32_t event{};
  std::uint32_t telescope{};
  bool weighted{};
  // Spherical atmosphere (atmospheric_acceptance_curved) if set, plane-parallel otherwise.
  std::optional<CurvedPath> curved{};
  // Sphere rows (frame G, sphere representation): atmospheric_acceptance_sphere, which takes
  // precedence over curved. The chunk must then carry the cz and path columns.
  std::optional<SphereGeometry> sphere{};
  // Threads for the rows of the chunk (1: in place). The result does not depend on it.
  unsigned threads{1};
};

namespace detail {

// Validation shared by the ground-frame steps: frame G, direction cosines, bunch weights,
// optional cz, zem (when needed), wavelength values. Weights are at most max_bunch_photons
// only when they are thinned (limit_weights); weighted output keeps any finite weight >= 0.
// Throws std::invalid_argument.
inline void validate_ground_chunk(obdeect::stream::PhotonChunk &chunk, bool need_zem,
                                  bool limit_weights) {
  using obdeect::stream::Column;
  if (chunk.frame != obdeect::stream::Frame::ground) {
    throw std::invalid_argument("photons must be in frame ground (G)");
  }
  chunk.validate();
  const std::size_t rows = chunk.rows();
  if (rows == 0) {
    return;
  }
  if (need_zem && !chunk.has(Column::zem)) {
    throw std::invalid_argument("extinction needs the zem column");
  }
  if (chunk.has(Column::p)) {
    throw std::invalid_argument("input already has an acceptance column p");
  }
  if (!chunk.has(Column::wavelength)) {
    chunk[Column::wavelength].assign(rows, 0.0F);
  }

  const auto &cx = chunk[Column::cx];
  const auto &cy = chunk[Column::cy];
  const auto &w = chunk[Column::w];
  for (std::size_t row = 0; row < rows; ++row) {
    if (!obdeect::stream::direction_in_unit_disc(cx[row], cy[row])) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    if (chunk.has(Column::cz)) {
      const auto cz_given = static_cast<double>(chunk[Column::cz][row]);
      // cz of either sign (upward photons, docs/propagation/propagation.md). cz_is_consistent is
      // false for NaN; the finiteness test is kept as the first check.
      if (!std::isfinite(cz_given) ||
          !obdeect::stream::cz_is_consistent(cz_given, cx[row], cy[row])) {
        throw std::invalid_argument(
            "cz column must be finite and match cx, cy (cz = +/-sqrt(1-cx^2-cy^2))");
      }
    }
    if (!(w[row] >= 0.0F && std::isfinite(w[row]))) {
      throw std::invalid_argument("bunch weight must be finite and >= 0");
    }
    if (limit_weights && static_cast<double>(w[row]) > static_cast<double>(max_bunch_photons)) {
      throw std::invalid_argument("bunch weight above max_bunch_photons (2048) cannot be "
                                  "thinned: use weighted mode or a smaller bunch size");
    }
    if (need_zem && !std::isfinite(chunk[Column::zem][row])) {
      throw std::invalid_argument("zem must be finite");
    }
  }

  const auto &wavelength = chunk[Column::wavelength];
  for (std::size_t row = 0; row < rows; ++row) {
    // 0 means not assigned; any other value must be finite and non-negative.
    if (wavelength[row] != 0.0F) {
      check_wavelength(static_cast<double>(wavelength[row]));
    }
  }
}

} // namespace detail

// Geometry-only step: the same checks and wavelength assignment as apply_extinction_chunk, but
// no extinction, envelope or thinning. Every row is kept with its weight; zem is not needed.
[[nodiscard]] inline std::vector<obdeect::stream::PhotonChunk>
assign_wavelength_chunk(obdeect::stream::PhotonChunk chunk, const WavelengthRange &wavelengths,
                        std::uint64_t seed, std::uint32_t event, std::uint32_t telescope) {
  wavelengths.validate();
  detail::validate_ground_chunk(chunk, false, false);
  if (chunk.rows() == 0) {
    return {};
  }
  assign_wavelengths(chunk, wavelengths, seed, event, telescope);
  std::vector<obdeect::stream::PhotonChunk> out;
  out.push_back(std::move(chunk));
  return out;
}

// Process one chunk (at most 2^16 rows). first_index of the chunk is the index of its
// first photon among all photons of this telescope and event; out_first_index numbers the
// output photons (thinned mode). Returns the output chunks (thinned mode may return more
// than one chunk; weighted mode returns the input chunk).
[[nodiscard]] inline std::vector<obdeect::stream::PhotonChunk>
apply_extinction_chunk(obdeect::stream::PhotonChunk chunk, const ExtinctionSettings &settings,
                       std::uint64_t out_first_index) {
  using obdeect::stream::Column;
  if (settings.table == nullptr || settings.depths == nullptr || settings.envelope == nullptr) {
    throw std::invalid_argument("extinction settings incomplete");
  }
  settings.wavelengths.validate();
  detail::validate_ground_chunk(chunk, true, !settings.weighted);
  const std::size_t rows = chunk.rows();
  if (rows == 0) {
    return {};
  }
  // Each stage runs over all rows before the next one, so the first error is the same for any
  // thread count. Rows are independent within a stage.
  const unsigned threads = settings.threads;
  const auto &wavelength = chunk[Column::wavelength];
  for_row_blocks(rows, threads, [&](std::size_t begin, std::size_t end) {
    assign_wavelengths(chunk, settings.wavelengths, settings.seed, settings.event,
                       settings.telescope, begin, end);
  });

  // Acceptance p = p_atm * E(lambda); rows at or above MAX_LAMBDA are dropped (p = 0).
  std::vector<double> p(rows, 0.0);
  for_row_blocks(rows, threads, [&](std::size_t begin, std::size_t end) {
    if (settings.sphere) {
      atmospheric_acceptance_sphere<double>(chunk, *settings.table, *settings.depths,
                                            *settings.sphere, p, begin, end);
    } else if (settings.curved) {
      atmospheric_acceptance_curved<double>(chunk, *settings.table, *settings.depths,
                                            *settings.curved, p, begin, end);
    } else {
      atmospheric_acceptance<double>(chunk, *settings.table, *settings.depths, p, begin, end);
    }
  });
  const Envelope &envelope = *settings.envelope;
  const bool linear = settings.table->wavelength_interpolation() == WavelengthInterpolation::linear;
  for_row_blocks(rows, threads, [&](std::size_t begin, std::size_t end) {
    for (std::size_t row = begin; row < end; ++row) {
      const int iwl = table_wavelength(static_cast<double>(wavelength[row]));
      if (iwl < 0) {
        p[row] = 0.0;
        continue;
      }
      // Floor mode: E at the integer nm of the row. Linear mode: E at the continuous wavelength,
      // linear between the nodes of the envelope (as the table lookups).
      const double lambda = linear ? static_cast<double>(wavelength[row]) : iwl;
      if (lambda < envelope.min_wavelength_nm() || lambda > envelope.max_wavelength_nm()) {
        throw std::invalid_argument("wavelength outside the envelope table");
      }
      p[row] *= envelope(lambda);
    }
  });

  if (settings.weighted) {
    weight_chunk(chunk, p, true);
    std::vector<obdeect::stream::PhotonChunk> out;
    out.push_back(std::move(chunk));
    return out;
  }
  return thin_chunk(chunk, p, settings.seed, settings.event, settings.telescope, out_first_index,
                    true, stream_thinning, threads);
}

// Envelope step for scattered photons (docs/scattering.md, frame G). Their transmission from the
// scattering point to the telescope sphere is in w, so the acceptance is E(lambda) times
// exp(-od_focus(lambda)), the in-telescope part of the extinction (the same od_focus, and the same
// wavelength mode, as apply_extinction_chunk). Thinned or weighted as apply_extinction_chunk;
// wavelengths must be assigned.
[[nodiscard]] inline std::vector<obdeect::stream::PhotonChunk>
apply_envelope_chunk(obdeect::stream::PhotonChunk chunk, const Envelope &envelope,
                     const TransmissionTable &table, const TelescopeDepths &depths, bool weighted,
                     std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                     std::uint64_t out_first_index) {
  using obdeect::stream::Column;
  if (chunk.frame != obdeect::stream::Frame::ground) {
    throw std::invalid_argument("envelope step needs photons in frame ground (G)");
  }
  chunk.validate();
  const std::size_t rows = chunk.rows();
  if (rows == 0) {
    return {};
  }
  if (!chunk.has(Column::wavelength) || chunk.has(Column::p)) {
    throw std::invalid_argument("envelope step needs wavelengths and no acceptance column p");
  }
  const auto &wavelength = chunk[Column::wavelength];
  const bool linear = table.wavelength_interpolation() == WavelengthInterpolation::linear;
  std::vector<double> p(rows, 0.0);
  for (std::size_t row = 0; row < rows; ++row) {
    if (wavelength[row] == 0.0F) {
      throw std::invalid_argument("envelope step needs assigned wavelengths");
    }
    const int iwl = table_wavelength(static_cast<double>(wavelength[row]));
    if (iwl < 0) {
      continue; // at or above MAX_LAMBDA: p = 0
    }
    const double lambda = linear ? static_cast<double>(wavelength[row]) : iwl;
    if (lambda < envelope.min_wavelength_nm() || lambda > envelope.max_wavelength_nm()) {
      throw std::invalid_argument("wavelength outside the envelope table");
    }
    // od_focus is 0 for a path of 0 m, so the factor is exactly 1 then.
    const double od_focus =
        table.telescope_depths_at(depths, static_cast<double>(wavelength[row])).od_focus;
    p[row] = envelope(lambda) * std::exp(-od_focus);
  }
  if (weighted) {
    weight_chunk(chunk, p, true);
    std::vector<obdeect::stream::PhotonChunk> out;
    out.push_back(std::move(chunk));
    return out;
  }
  return thin_chunk(chunk, p, seed, event, telescope, out_first_index, true,
                    stream_scattered_thinning);
}

} // namespace obdeect::propagation
