#pragma once

// One chunk of the atmosphere step for telescope photons in frame ground:
// wavelength assignment, atmospheric acceptance, envelope, then weighted or thinned output.

#include "dlufd/extinction.hpp"
#include "dlufd/streams.hpp"
#include "dlufd/thinning.hpp"
#include "dlufd/transmission.hpp"
#include "dlufd/wavelength.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace dlufd {

struct ExtinctionSettings {
  const TransmissionTable *table{};
  const TelescopeDepths *depths{};
  const Envelope *envelope{};
  WavelengthRange wavelengths{};
  std::uint64_t seed{};
  std::uint32_t event{};
  std::uint32_t telescope{};
  bool weighted{};
};

// Process one chunk (at most 2^16 rows). first_index of the chunk is the index of its
// first photon among all photons of this telescope and event; out_first_index numbers the
// output photons (thinned mode). Returns the output chunks (thinned mode may return more
// than one chunk; weighted mode returns the input chunk).
[[nodiscard]] inline std::vector<gstrickd::PhotonChunk>
apply_extinction_chunk(gstrickd::PhotonChunk chunk, const ExtinctionSettings &settings,
                       std::uint64_t out_first_index) {
  using gstrickd::Column;
  if (settings.table == nullptr || settings.depths == nullptr || settings.envelope == nullptr) {
    throw std::invalid_argument("extinction settings incomplete");
  }
  settings.wavelengths.validate();
  if (chunk.frame != gstrickd::Frame::ground) {
    throw std::invalid_argument("extinction needs photons in frame ground (G)");
  }
  chunk.validate();
  const std::size_t rows = chunk.rows();
  if (rows == 0) {
    return {};
  }
  if (!chunk.has(Column::zem)) {
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
  const auto &zem = chunk[Column::zem];
  for (std::size_t row = 0; row < rows; ++row) {
    const double r2 =
        static_cast<double>(cx[row]) * cx[row] + static_cast<double>(cy[row]) * cy[row];
    if (!(r2 <= 1.0)) {
      throw std::invalid_argument("direction cosines must satisfy cx^2 + cy^2 <= 1");
    }
    if (chunk.has(Column::cz)) {
      const double cz = -std::sqrt(1.0 - r2);
      const auto cz_given = static_cast<double>(chunk[Column::cz][row]);
      // NaN fails the finiteness test; a NaN difference would pass the tolerance test.
      if (!std::isfinite(cz_given) || std::abs(cz_given - cz) > 1e-5) {
        throw std::invalid_argument(
            "cz column must be finite and match cx, cy (cz = -sqrt(1-cx^2-cy^2))");
      }
    }
    if (!(w[row] >= 0.0F &&
          static_cast<double>(w[row]) <= static_cast<double>(max_bunch_photons))) {
      throw std::invalid_argument("bunch weight must be within [0, max_bunch_photons]");
    }
    if (!std::isfinite(zem[row])) {
      throw std::invalid_argument("zem must be finite");
    }
  }

  auto &wavelength = chunk[Column::wavelength];
  for (std::size_t row = 0; row < rows; ++row) {
    // 0 means not assigned; any other value must be a finite, non-negative integer.
    if (wavelength[row] != 0.0F) {
      check_wavelength(static_cast<double>(wavelength[row]));
    }
  }
  assign_wavelengths(chunk, settings.wavelengths, settings.seed, settings.event,
                     settings.telescope);

  // Acceptance p = p_atm * E(lambda); rows at or above MAX_LAMBDA are dropped (p = 0).
  std::vector<double> p(rows, 0.0);
  atmospheric_acceptance<double>(chunk, *settings.table, *settings.depths, p);
  const Envelope &envelope = *settings.envelope;
  for (std::size_t row = 0; row < rows; ++row) {
    const int iwl = table_wavelength(static_cast<double>(wavelength[row]));
    if (iwl < 0) {
      p[row] = 0.0;
      continue;
    }
    if (iwl < envelope.min_wavelength_nm() || iwl > envelope.max_wavelength_nm()) {
      throw std::invalid_argument("wavelength outside the envelope table");
    }
    p[row] *= envelope(static_cast<double>(iwl));
  }

  if (settings.weighted) {
    weight_chunk(chunk, p, true);
    std::vector<gstrickd::PhotonChunk> out;
    out.push_back(std::move(chunk));
    return out;
  }
  return thin_chunk(chunk, p, settings.seed, settings.event, settings.telescope, out_first_index,
                    true);
}

// Envelope step for scattered photons (docs/scattering.md, frame G). Their transmission to the
// telescope is in w, so the acceptance is the envelope E(lambda) only (no atmospheric
// extinction). Thinned or weighted as apply_extinction_chunk; wavelengths must be assigned.
[[nodiscard]] inline std::vector<gstrickd::PhotonChunk>
apply_envelope_chunk(gstrickd::PhotonChunk chunk, const Envelope &envelope, bool weighted,
                     std::uint64_t seed, std::uint32_t event, std::uint32_t telescope,
                     std::uint64_t out_first_index) {
  using gstrickd::Column;
  if (chunk.frame != gstrickd::Frame::ground) {
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
  std::vector<double> p(rows, 0.0);
  for (std::size_t row = 0; row < rows; ++row) {
    if (wavelength[row] == 0.0F) {
      throw std::invalid_argument("envelope step needs assigned wavelengths");
    }
    const int iwl = table_wavelength(static_cast<double>(wavelength[row]));
    if (iwl < 0) {
      continue; // at or above MAX_LAMBDA: p = 0
    }
    if (iwl < envelope.min_wavelength_nm() || iwl > envelope.max_wavelength_nm()) {
      throw std::invalid_argument("wavelength outside the envelope table");
    }
    p[row] = envelope(static_cast<double>(iwl));
  }
  if (weighted) {
    weight_chunk(chunk, p, true);
    std::vector<gstrickd::PhotonChunk> out;
    out.push_back(std::move(chunk));
    return out;
  }
  return thin_chunk(chunk, p, seed, event, telescope, out_first_index, true,
                    stream_scattered_thinning);
}

} // namespace dlufd
