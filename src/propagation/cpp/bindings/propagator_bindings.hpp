#pragma once

// Python API of the atmosphere stage: dlufd.propagate on stream files (the stream I/O is
// the C++ core's). Included by core.cpp after the nanobind headers.

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include "dlufd/propagator.hpp"

#include <cstdint>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace nb = nanobind;

inline void bind_propagator(nb::module_ &module) {
  module.def(
      "propagate",
      [](const std::string &input_path, const std::string &output_path,
         const std::string &transmission_path, double site_altitude_m,
         const std::optional<std::string> &profile_path,
         const std::vector<std::pair<std::string, std::string>> &envelope_paths,
         const std::optional<std::tuple<double, double>> &wavelength_range, std::uint64_t seed,
         bool weighted, bool store_p, const std::optional<std::string> &scattering_path,
         std::uint32_t scattering_points) {
        dlufd::PropagatorConfig config;
        config.transmission.emplace(
            dlufd::TransmissionTable::from_file(transmission_path, site_altitude_m));
        config.transmission_source = transmission_path;
        config.site_altitude_m = site_altitude_m;
        if (profile_path) {
          config.profile.emplace(dlufd::AtmosphereProfile::from_file(*profile_path));
          config.profile_source = *profile_path;
        }
        for (const auto &[name, path] : envelope_paths) {
          config.envelope.push_back(dlufd::read_envelope_file(name, path));
        }
        if (wavelength_range) {
          config.wavelength_range = dlufd::WavelengthRange{std::get<0>(*wavelength_range),
                                                           std::get<1>(*wavelength_range)};
        }
        config.seed = seed;
        config.mode = weighted ? dlufd::OutputMode::weighted : dlufd::OutputMode::thinned;
        config.store_p = store_p;
        if (scattering_path) {
          config.scattering.enabled = true;
          config.scattering.points_per_photon = scattering_points;
          config.scattering_table.emplace(dlufd::ScatteringTable::from_file(*scattering_path));
          config.scattering_source = *scattering_path;
        }

        std::ifstream in(input_path, std::ios::binary);
        if (!in) {
          throw std::runtime_error("cannot open input stream: " + input_path);
        }
        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
          throw std::runtime_error("cannot open output stream: " + output_path);
        }
        const dlufd::PropagationSummary summary =
            dlufd::propagate_stream(in, out, std::move(config));
        out.close();
        nb::dict result;
        result["photons_in"] = summary.photons_in;
        result["photons_out"] = summary.photons_out;
        result["dropped"] = summary.dropped;
        result["weight_out"] = summary.weight_out;
        return result;
      },
      nb::arg("input_path"), nb::arg("output_path"), nb::arg("transmission_path"),
      nb::arg("site_altitude_m"), nb::arg("profile_path") = nb::none(),
      nb::arg("envelope_paths") = std::vector<std::pair<std::string, std::string>>{},
      nb::arg("wavelength_range") = nb::none(), nb::arg("seed"), nb::arg("weighted") = false,
      nb::arg("store_p") = false, nb::arg("scattering_path") = nb::none(),
      nb::arg("scattering_points") = std::uint32_t{8},
      "Run the atmosphere stage on a photon stream file. envelope_paths: (name, file) pairs of "
      "two-column text tables (nm, value). scattering_path: optional scattering table (single "
      "scattering of the primaries, off when None).");
}
