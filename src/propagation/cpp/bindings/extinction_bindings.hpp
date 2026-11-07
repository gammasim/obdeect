#pragma once

// nanobind bindings of the atmosphere step: TransmissionTable and apply_extinction.
// Included by core.cpp only; the C++ core stays free of Python.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include "dlufd/extinction.hpp"
#include "dlufd/pipeline.hpp"
#include "dlufd/thinning.hpp"
#include "dlufd/transmission.hpp"
#include "dlufd/wavelength.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace dlufd_bindings {

namespace nb = nanobind;

using DoubleArray = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using FloatArray = nb::ndarray<const float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// Hand a vector to Python without copying; the returned array owns the vector.
template <class T> nb::ndarray<nb::numpy, T, nb::ndim<1>> to_numpy(std::vector<T> values) {
  if (values.capacity() == 0) {
    values.reserve(1); // a non-null data pointer also for empty arrays
  }
  auto *holder = new std::vector<T>(std::move(values));
  nb::capsule owner(holder, [](void *p) noexcept { delete static_cast<std::vector<T> *>(p); });
  return nb::ndarray<nb::numpy, T, nb::ndim<1>>(holder->data(), {holder->size()}, owner);
}

inline std::vector<double> to_vector(const DoubleArray &array) {
  return std::vector<double>(array.data(), array.data() + array.size());
}

inline void bind_table(nb::module_ &module) {
  using dlufd::TransmissionTable;
  nb::class_<TransmissionTable>(module, "TransmissionTable")
      .def(nb::new_([](const std::string &path, double site_altitude_m) {
             return TransmissionTable::from_file(path, site_altitude_m);
           }),
           nb::arg("path"), nb::arg("site_altitude_m"))
      .def_prop_ro("observation_level_m", &TransmissionTable::observation_level_m)
      .def_prop_ro("first_wavelength_nm", &TransmissionTable::first_wavelength_nm)
      .def_prop_ro("last_wavelength_nm", &TransmissionTable::last_wavelength_nm)
      .def(
          "optical_depth",
          [](const TransmissionTable &self, const DoubleArray &wavelength_nm,
             const DoubleArray &height_m) {
            if (height_m.size() != wavelength_nm.size()) {
              throw std::invalid_argument("arrays must have the same length");
            }
            std::vector<double> out(wavelength_nm.size());
            for (std::size_t i = 0; i < out.size(); ++i) {
              out[i] = self.optical_depth<double>(wavelength_nm.data()[i], height_m.data()[i]);
            }
            return to_numpy(std::move(out));
          },
          nb::arg("wavelength_nm"), nb::arg("height_m"))
      .def(
          "transmission",
          [](const TransmissionTable &self, const DoubleArray &wavelength_nm,
             const DoubleArray &zem_m, const DoubleArray &airmass) {
            const std::size_t n = wavelength_nm.size();
            if (zem_m.size() != n || airmass.size() != n) {
              throw std::invalid_argument("arrays must have the same length");
            }
            std::vector<double> out(n);
            for (std::size_t i = 0; i < n; ++i) {
              out[i] = self.transmission<double>(wavelength_nm.data()[i], zem_m.data()[i],
                                                 airmass.data()[i]);
            }
            return to_numpy(std::move(out));
          },
          nb::arg("wavelength_nm"), nb::arg("zem_m"), nb::arg("airmass"))
      .def(
          "transmission2",
          [](const TransmissionTable &self, const DoubleArray &wavelength_nm,
             const DoubleArray &zem_m, const DoubleArray &airmass, const DoubleArray &od_tel,
             const DoubleArray &od_focus) {
            const std::size_t n = wavelength_nm.size();
            if (zem_m.size() != n || airmass.size() != n || od_tel.size() != n ||
                od_focus.size() != n) {
              throw std::invalid_argument("arrays must have the same length");
            }
            std::vector<double> out(n);
            for (std::size_t i = 0; i < n; ++i) {
              out[i] = self.transmission2<double>(wavelength_nm.data()[i], zem_m.data()[i],
                                                  airmass.data()[i], od_tel.data()[i],
                                                  od_focus.data()[i]);
            }
            return to_numpy(std::move(out));
          },
          nb::arg("wavelength_nm"), nb::arg("zem_m"), nb::arg("airmass"), nb::arg("od_tel"),
          nb::arg("od_focus"))
      .def(
          "telescope_depths",
          [](const TransmissionTable &self, double telescope_altitude_m, double pathlen_m) {
            const auto depths = self.telescope_depths(telescope_altitude_m, pathlen_m);
            return nb::make_tuple(
                to_numpy(std::vector<double>(depths.od_tel.begin(), depths.od_tel.end())),
                to_numpy(std::vector<double>(depths.od_focus.begin(), depths.od_focus.end())));
          },
          nb::arg("telescope_altitude_m"), nb::arg("pathlen_m") = 0.0);
}

// Map a column name of the Python API to its gstrickd column.
inline gstrickd::Column column_from_name(const std::string &name) {
  for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
    const auto column = static_cast<gstrickd::Column>(c);
    if (name == gstrickd::column_name(column)) {
      return column;
    }
  }
  throw std::invalid_argument("unknown photon column: " + name);
}

// Columns given as equal-length float32 arrays, keyed by column name.
inline gstrickd::PhotonChunk chunk_from_columns(const nb::dict &columns,
                                                std::uint64_t first_index) {
  std::vector<std::pair<gstrickd::Column, FloatArray>> arrays;
  std::size_t rows = 0;
  std::uint32_t mask = 0;
  for (auto item : columns) {
    const auto column = column_from_name(nb::cast<std::string>(item.first));
    const auto array = nb::cast<FloatArray>(item.second);
    if (!arrays.empty() && array.size() != rows) {
      throw std::invalid_argument("photon columns differ in length");
    }
    rows = array.size();
    mask |= gstrickd::column_bit(column);
    arrays.emplace_back(column, array);
  }
  gstrickd::PhotonChunk chunk;
  chunk.first_index = first_index;
  chunk.resize(rows, mask);
  for (const auto &[column, array] : arrays) {
    std::copy(array.data(), array.data() + array.size(), chunk[column].begin());
  }
  return chunk;
}

inline void bind_apply(nb::module_ &module) {
  module.def(
      "apply_extinction",
      [](const nb::dict &columns, const dlufd::TransmissionTable &table,
         double telescope_altitude_m, const DoubleArray &envelope_wavelength,
         const DoubleArray &envelope_value, std::uint64_t seed, std::uint32_t event,
         std::uint32_t telescope, std::uint64_t first_index, bool weighted,
         std::tuple<double, double> wavelength_range, double pathlen_m) {
        const dlufd::Envelope envelope(to_vector(envelope_wavelength), to_vector(envelope_value));
        const auto depths = table.telescope_depths(telescope_altitude_m, pathlen_m);
        dlufd::ExtinctionSettings settings;
        settings.table = &table;
        settings.depths = &depths;
        settings.envelope = &envelope;
        settings.wavelengths = {std::get<0>(wavelength_range), std::get<1>(wavelength_range)};
        settings.seed = seed;
        settings.event = event;
        settings.telescope = telescope;
        settings.weighted = weighted;

        // The input is processed in slices of at most 2^16 rows. The photon index and the
        // output index continue across slices, so the result does not depend on slicing.
        const gstrickd::PhotonChunk input = chunk_from_columns(columns, first_index);
        std::vector<std::vector<float>> out(gstrickd::n_columns);
        std::uint64_t out_index = 0;
        for (std::size_t offset = 0; offset < input.rows(); offset += gstrickd::max_chunk_rows) {
          const std::size_t rows =
              std::min<std::size_t>(gstrickd::max_chunk_rows, input.rows() - offset);
          gstrickd::PhotonChunk slice;
          slice.telescope_id = input.telescope_id;
          slice.frame = input.frame;
          slice.first_index = first_index + offset;
          slice.resize(rows, input.mask());
          for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
            if (input.has(static_cast<gstrickd::Column>(c))) {
              const auto begin = input.columns[c].begin() + static_cast<std::ptrdiff_t>(offset);
              std::copy(begin, begin + static_cast<std::ptrdiff_t>(rows), slice.columns[c].begin());
            }
          }
          for (auto &chunk : dlufd::apply_extinction_chunk(std::move(slice), settings, out_index)) {
            for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
              if (chunk.has(static_cast<gstrickd::Column>(c))) {
                out[c].insert(out[c].end(), chunk.columns[c].begin(), chunk.columns[c].end());
              }
            }
            out_index += chunk.rows();
          }
        }
        nb::dict result;
        for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
          const auto column = static_cast<gstrickd::Column>(c);
          if (!out[c].empty()) {
            result[gstrickd::column_name(column)] = nb::cast(to_numpy(std::move(out[c])));
          }
        }
        return result;
      },
      nb::arg("columns"), nb::arg("table"), nb::arg("telescope_altitude_m"),
      nb::arg("envelope_wavelength"), nb::arg("envelope_value"), nb::arg("seed"), nb::arg("event"),
      nb::arg("telescope"), nb::arg("first_index"), nb::arg("weighted"),
      nb::arg("wavelength_range"), nb::arg("pathlen_m") = 0.0);
}

} // namespace dlufd_bindings

// Register TransmissionTable and apply_extinction on the extension module.
inline void bind_extinction(nanobind::module_ &module) {
  dlufd_bindings::bind_table(module);
  dlufd_bindings::bind_apply(module);
}
