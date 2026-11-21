#pragma once

// Python API of `obdeect propagation to-sky`: obdeect.propagation.to_sky on stream files, the
// same conversion as the command-line tool. Included by core.cpp after the nanobind headers.

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>

#include "obdeect/propagation/ground_to_sky.hpp"
#include "obdeect/propagation/profile.hpp"
#include "obdeect/propagation/propagator.hpp"

#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>

namespace nb = nanobind;

inline void bind_to_sky(nb::module_ &module) {
  module.def(
      "to_sky",
      [](const std::string &input_path, const std::string &output_path,
         const std::optional<std::string> &profile_path) {
        std::optional<obdeect::propagation::AtmosphereProfile> profile;
        if (profile_path) {
          profile.emplace(obdeect::propagation::AtmosphereProfile::from_file(*profile_path));
        }
        if (obdeect::propagation::same_file(input_path, output_path)) {
          throw std::invalid_argument("input and output are the same file: " + input_path);
        }
        std::ifstream in(input_path, std::ios::binary);
        if (!in) {
          throw std::runtime_error("cannot open input stream: " + input_path);
        }
        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
          throw std::runtime_error("cannot open output stream: " + output_path);
        }
        const obdeect::propagation::ToSkySummary summary =
            obdeect::propagation::ground_to_sky_stream(in, out, profile,
                                                       "obdeect propagation to-sky");
        out.close();
        nb::dict result;
        result["photons_in"] = summary.photons_in;
        result["photons_out"] = summary.photons_out;
        result["dropped"] = summary.dropped;
        return result;
      },
      nb::arg("input_path"), nb::arg("output_path"), nb::arg("profile_path") = nb::none(),
      "Convert the frame G photons of a stream to frame sky (the inverse of the atmosphere stage's "
      "sky-to-G step; docs/propagation/ground-to-sky.md). profile_path: the run's atmospheric "
      "profile (the air delay is inverted; the run needs refraction_applied = 1 with a profile). "
      "Returns photons_in, photons_out and dropped.");
}
