#pragma once

// nanobind bindings of single scattering: ScatteringTable and scatter.
// Included by core.cpp only; the C++ core stays free of Python.

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include "extinction_bindings.hpp"

#include "dlufd/scattering.hpp"
#include "dlufd/transmission.hpp"

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

using Telescope = std::tuple<std::uint32_t, double, double, double, double>;

inline void bind_scattering_table(nb::module_ &module) {
  using dlufd::ScatteringTable;
  nb::class_<ScatteringTable>(module, "ScatteringTable")
      .def(nb::new_([](const std::string &path) { return ScatteringTable::from_file(path); }),
           nb::arg("path"))
      .def_prop_ro("first_wavelength_nm", &ScatteringTable::first_wavelength_nm)
      .def_prop_ro("last_wavelength_nm", &ScatteringTable::last_wavelength_nm)
      .def_prop_ro("min_altitude_m", &ScatteringTable::min_altitude_m)
      .def_prop_ro("max_altitude_m", &ScatteringTable::max_altitude_m)
      .def(
          "profile",
          [](const ScatteringTable &self, const DoubleArray &wavelength_nm,
             const DoubleArray &altitude_m) {
            if (altitude_m.size() != wavelength_nm.size()) {
              throw std::invalid_argument("arrays must have the same length");
            }
            const std::size_t n = wavelength_nm.size();
            std::vector<double> tau_mol(n), tau_aer(n), ssa(n), asymmetry(n);
            for (std::size_t i = 0; i < n; ++i) {
              const auto s = self.state<double>(wavelength_nm.data()[i], altitude_m.data()[i]);
              tau_mol[i] = s.tau_mol;
              tau_aer[i] = s.tau_aer;
              ssa[i] = s.ssa;
              asymmetry[i] = s.asymmetry;
            }
            return nb::make_tuple(to_numpy(std::move(tau_mol)), to_numpy(std::move(tau_aer)),
                                  to_numpy(std::move(ssa)), to_numpy(std::move(asymmetry)));
          },
          nb::arg("wavelength_nm"), nb::arg("altitude_m"));
}

inline void bind_scatter(nb::module_ &module) {
  module.def(
      "scatter",
      [](const nb::dict &columns, const dlufd::ScatteringTable &table,
         const dlufd::TransmissionTable &transmission, const std::vector<Telescope> &telescopes,
         std::uint64_t seed, std::uint32_t event, std::uint64_t first_index,
         std::uint32_t n_points) {
        std::vector<dlufd::ScatterTelescope> sphere;
        for (const auto &[id, x, y, z, radius] : telescopes) {
          sphere.push_back({id, x, y, z, radius});
        }
        dlufd::ScatteringConfig config;
        config.enabled = true;
        config.points_per_photon = n_points;

        // Input in slices of at most 2^16 rows; the photon index continues across slices and
        // the output of each slice is appended in order, so the result does not depend on it.
        gstrickd::PhotonChunk input = chunk_from_columns(columns, first_index);
        input.frame = gstrickd::Frame::sky; // primary photons of the sky frame
        std::vector<std::vector<gstrickd::PhotonChunk>> lists(sphere.size());
        for (std::size_t offset = 0; offset < input.rows(); offset += gstrickd::max_chunk_rows) {
          const std::size_t rows =
              std::min<std::size_t>(gstrickd::max_chunk_rows, input.rows() - offset);
          gstrickd::PhotonChunk slice;
          slice.frame = input.frame;
          slice.first_index = first_index + offset;
          slice.resize(rows, input.mask());
          for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
            if (input.has(static_cast<gstrickd::Column>(c))) {
              const auto begin = input.columns[c].begin() + static_cast<std::ptrdiff_t>(offset);
              std::copy(begin, begin + static_cast<std::ptrdiff_t>(rows), slice.columns[c].begin());
            }
          }
          dlufd::scatter_chunk(config, table, transmission, sphere, slice, seed, event,
                               first_index + offset, lists);
        }

        // One dict of float32 columns per telescope, keyed by telescope id.
        nb::dict result;
        for (std::size_t i = 0; i < sphere.size(); ++i) {
          std::vector<std::vector<float>> out(gstrickd::n_columns);
          for (const auto &chunk : lists[i]) {
            for (const auto column : dlufd::scatter_output_columns) {
              const auto c = static_cast<std::uint32_t>(column);
              out[c].insert(out[c].end(), chunk.columns[c].begin(), chunk.columns[c].end());
            }
          }
          nb::dict photons;
          for (const auto column : dlufd::scatter_output_columns) {
            const auto c = static_cast<std::uint32_t>(column);
            photons[gstrickd::column_name(column)] = nb::cast(to_numpy(std::move(out[c])));
          }
          result[nb::int_(sphere[i].id)] = photons;
        }
        return result;
      },
      nb::arg("columns"), nb::arg("table"), nb::arg("transmission"), nb::arg("telescopes"),
      nb::arg("seed"), nb::arg("event"), nb::arg("first_index"), nb::arg("n_points"));
}

} // namespace dlufd_bindings

// Register ScatteringTable and scatter on the extension module.
inline void bind_scattering(nanobind::module_ &module) {
  dlufd_bindings::bind_scattering_table(module);
  dlufd_bindings::bind_scatter(module);
}
