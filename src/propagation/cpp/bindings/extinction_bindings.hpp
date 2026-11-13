#pragma once

// nanobind bindings of the atmosphere step: TransmissionTable and apply_extinction.
// Included by core.cpp only; the C++ core stays free of Python.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include "obdeect/propagation/extinction.hpp"
#include "obdeect/propagation/pipeline.hpp"
#include "obdeect/propagation/thinning.hpp"
#include "obdeect/propagation/transmission.hpp"
#include "obdeect/propagation/wavelength.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace obdeect_propagation_bindings {

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
  using obdeect::propagation::TransmissionTable;
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
inline obdeect::stream::Column column_from_name(const std::string &name) {
  for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c) {
    const auto column = static_cast<obdeect::stream::Column>(c);
    if (name == obdeect::stream::column_name(column)) {
      return column;
    }
  }
  throw std::invalid_argument("unknown photon column: " + name);
}

// Input columns given as equal-length float32 arrays, keyed by column name. The arrays
// are not copied; mask is the set of columns given (also when there are no rows).
struct InputColumns {
  std::vector<std::pair<obdeect::stream::Column, FloatArray>> arrays;
  std::size_t rows{};
  std::uint32_t mask{};
};

inline InputColumns input_columns(const nb::dict &columns) {
  InputColumns input;
  for (auto item : columns) {
    const auto column = column_from_name(nb::cast<std::string>(item.first));
    const auto array = nb::cast<FloatArray>(item.second);
    if (!input.arrays.empty() && array.size() != input.rows) {
      throw std::invalid_argument("photon columns differ in length");
    }
    input.rows = array.size();
    input.mask |= obdeect::stream::column_bit(column);
    input.arrays.emplace_back(column, array);
  }
  return input;
}

// Chunk of rows [offset, offset + count) of the input; first_index is the index of row
// offset among all photons. count must be at most obdeect::stream::max_chunk_rows.
inline obdeect::stream::PhotonChunk slice_chunk(const InputColumns &input, std::size_t offset,
                                         std::size_t count, std::uint64_t first_index) {
  obdeect::stream::PhotonChunk chunk;
  chunk.first_index = first_index;
  chunk.resize(count, input.mask);
  for (const auto &[column, array] : input.arrays) {
    const float *begin = array.data() + offset;
    std::copy(begin, begin + count, chunk[column].begin());
  }
  return chunk;
}

inline void bind_apply(nb::module_ &module) {
  module.def(
      "apply_extinction",
      [](const nb::dict &columns, const obdeect::propagation::TransmissionTable &table,
         double telescope_altitude_m, const DoubleArray &envelope_wavelength,
         const DoubleArray &envelope_value, std::uint64_t seed, std::uint32_t event,
         std::uint32_t telescope, std::uint64_t first_index, bool weighted,
         std::tuple<double, double> wavelength_range, double pathlen_m) {
        const obdeect::propagation::Envelope envelope(to_vector(envelope_wavelength), to_vector(envelope_value));
        const auto depths = table.telescope_depths(telescope_altitude_m, pathlen_m);
        obdeect::propagation::ExtinctionSettings settings;
        settings.table = &table;
        settings.depths = &depths;
        settings.envelope = &envelope;
        settings.wavelengths = {std::get<0>(wavelength_range), std::get<1>(wavelength_range)};
        settings.seed = seed;
        settings.event = event;
        settings.telescope = telescope;
        settings.weighted = weighted;

        // The input is processed in slices of at most 2^16 rows (the PhotonChunk limit). Each
        // slice keeps its photon index (first_index + offset), so the random draws do not
        // depend on slicing; the output index continues across slices.
        const InputColumns input = input_columns(columns);
        // Output columns: the input columns plus the wavelength and p columns, which the
        // pipeline always produces, also when no photon survives.
        const std::uint32_t output_mask = input.mask |
                                          obdeect::stream::column_bit(obdeect::stream::Column::wavelength) |
                                          obdeect::stream::column_bit(obdeect::stream::Column::p);
        std::vector<std::vector<float>> out(obdeect::stream::n_columns);
        std::uint64_t out_index = 0;
        for (std::size_t offset = 0; offset < input.rows; offset += obdeect::stream::max_chunk_rows) {
          const std::size_t rows =
              std::min<std::size_t>(obdeect::stream::max_chunk_rows, input.rows - offset);
          obdeect::stream::PhotonChunk slice = slice_chunk(input, offset, rows, first_index + offset);
          for (auto &chunk : obdeect::propagation::apply_extinction_chunk(std::move(slice), settings, out_index)) {
            for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c) {
              if (chunk.has(static_cast<obdeect::stream::Column>(c))) {
                out[c].insert(out[c].end(), chunk.columns[c].begin(), chunk.columns[c].end());
              }
            }
            out_index += chunk.rows();
          }
        }
        nb::dict result;
        for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c) {
          const auto column = static_cast<obdeect::stream::Column>(c);
          if ((output_mask >> c) & 1U) {
            result[obdeect::stream::column_name(column)] = nb::cast(to_numpy(std::move(out[c])));
          }
        }
        return result;
      },
      nb::arg("columns"), nb::arg("table"), nb::arg("telescope_altitude_m"),
      nb::arg("envelope_wavelength"), nb::arg("envelope_value"), nb::arg("seed"), nb::arg("event"),
      nb::arg("telescope"), nb::arg("first_index"), nb::arg("weighted"),
      nb::arg("wavelength_range"), nb::arg("pathlen_m") = 0.0);
}

} // namespace obdeect_propagation_bindings

// Register TransmissionTable and apply_extinction on the extension module.
inline void bind_extinction(nanobind::module_ &module) {
  obdeect_propagation_bindings::bind_table(module);
  obdeect_propagation_bindings::bind_apply(module);
}
