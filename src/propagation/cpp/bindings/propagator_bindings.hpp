#pragma once

// Python API of the atmosphere stage: obdeect.propagation.propagate on stream files (the stream I/O
// is the C++ core's). Included by core.cpp after the nanobind headers.
//
// The Python wrapper passes every argument: the defaults live in Python only.

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include "obdeect/propagation/propagator.hpp"
#include "stream_file.hpp"

#include <cstdint>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace obdeect::propagation {

// The dropped photons of a summary by reason (keys as DropCounts; the sum is "dropped").
inline nb::dict dropped_by_reason_dict(const DropCounts &counts) {
  nb::dict out;
  out["moving_away"] = counts.moving_away;
  out["horizontal"] = counts.horizontal;
  out["below_ground"] = counts.below_ground;
  out["turns_over"] = counts.turns_over;
  out["not_bendable"] = counts.not_bendable;
  out["misses_level"] = counts.misses_level;
  out["misses_sphere"] = counts.misses_sphere;
  return out;
}

// "thinned" or "weighted" (the Python mode argument); anything else is rejected.
inline OutputMode parse_output_mode(const std::string &text) {
  if (text == "thinned") {
    return OutputMode::thinned;
  }
  if (text == "weighted") {
    return OutputMode::weighted;
  }
  throw std::invalid_argument("mode must be thinned or weighted: '" + text + "'");
}

// 1-based telescope ids of a selection; the Python side gives integers.
inline std::vector<std::uint32_t> telescope_ids(const std::vector<std::int64_t> &ids,
                                                const char *what) {
  std::vector<std::uint32_t> selection;
  selection.reserve(ids.size());
  for (const std::int64_t id : ids) {
    if (id < 1 || id > static_cast<std::int64_t>(obdeect::stream::max_telescope)) {
      throw std::invalid_argument(std::string(what) + ": telescope ids are 1-based");
    }
    selection.push_back(static_cast<std::uint32_t>(id));
  }
  return selection;
}

} // namespace obdeect::propagation

inline void bind_propagator(nb::module_ &module) {
  module.def(
      "propagate",
      [](const std::string &input_path, const std::string &output_path,
         const std::optional<std::string> &transmission_path, double site_altitude_m,
         const std::optional<std::string> &profile_path,
         const std::vector<std::pair<std::string, std::string>> &envelope_paths,
         const std::optional<std::tuple<double, double>> &wavelength_range, std::uint64_t seed,
         const std::string &mode, bool store_p, const std::optional<std::string> &scattering_path,
         std::uint32_t scattering_points, double pathlen_m,
         const std::map<std::int64_t, double> &telescope_pathlen_m, bool geometry_only,
         const std::string &curvature, const std::optional<double> &earth_radius_m,
         const std::optional<std::vector<std::int64_t>> &telescopes, std::int64_t threads,
         const std::optional<std::vector<std::int64_t>> &scattering_primaries,
         const std::string &wavelength_interpolation,
         const std::optional<double> &dispersion_reference_nm,
         const std::string &ground_representation) {
        if (threads < 1 || threads > 256) {
          throw std::invalid_argument("threads must be an integer from 1 to 256");
        }
        obdeect::propagation::PropagatorConfig config;
        config.threads = static_cast<unsigned>(threads);
        config.geometry_only = geometry_only;
        if (telescopes) {
          config.telescopes = obdeect::propagation::telescope_ids(*telescopes, "telescopes");
        }
        if (scattering_primaries) {
          config.scattering_primaries =
              obdeect::propagation::telescope_ids(*scattering_primaries, "scattering_primaries");
        }
        config.curvature = obdeect::propagation::parse_curvature(curvature);
        config.ground_representation =
            obdeect::propagation::parse_ground_representation(ground_representation);
        if (earth_radius_m) {
          config.earth_radius_m = *earth_radius_m;
        }
        config.default_pathlen_m = pathlen_m;
        for (const auto &[id, metres] : telescope_pathlen_m) {
          if (id < 1 || id > static_cast<std::int64_t>(obdeect::stream::max_telescope)) {
            throw std::invalid_argument("telescope_pathlen_m: telescope ids are 1-based");
          }
          config.telescope_pathlen_m[static_cast<std::uint32_t>(id)] = metres;
        }
        if (transmission_path) {
          config.transmission.emplace(obdeect::propagation::TransmissionTable::from_file(
              *transmission_path, site_altitude_m,
              obdeect::propagation::parse_wavelength_interpolation(wavelength_interpolation)));
          config.transmission_source = *transmission_path;
        }
        config.dispersion_reference_nm = dispersion_reference_nm;
        config.site_altitude_m = site_altitude_m;
        if (profile_path) {
          config.profile.emplace(obdeect::propagation::AtmosphereProfile::from_file(*profile_path));
          config.profile_source = *profile_path;
        }
        for (const auto &[name, path] : envelope_paths) {
          config.envelope.push_back(obdeect::propagation::read_envelope_file(name, path));
        }
        if (wavelength_range) {
          config.wavelength_range = obdeect::propagation::WavelengthRange{
              std::get<0>(*wavelength_range), std::get<1>(*wavelength_range)};
        }
        config.seed = seed;
        config.mode = obdeect::propagation::parse_output_mode(mode);
        config.store_p = store_p;
        if (scattering_path) {
          config.scattering.points_per_photon = scattering_points;
          config.scattering_table.emplace(
              obdeect::propagation::ScatteringTable::from_file(*scattering_path));
          config.scattering_source = *scattering_path;
        }

        const obdeect::propagation::PropagationSummary summary =
            obdeect_propagation_bindings::convert_stream_file(
                input_path, output_path, [&config](std::istream &in, std::ostream &out) {
                  return obdeect::propagation::propagate_stream(in, out, std::move(config));
                });
        nb::dict result;
        result["photons_in"] = summary.photons_in;
        result["photons_out"] = summary.photons_out;
        result["dropped"] = summary.dropped;
        result["dropped_by_reason"] =
            obdeect::propagation::dropped_by_reason_dict(summary.dropped_by_reason);
        result["weight_out"] = summary.weight_out;
        result["thinned_out"] = summary.thinned_out;
        result["scattered_dropped"] = summary.scattered_dropped;
        result["scattered_dropped_by_reason"] =
            obdeect::propagation::dropped_by_reason_dict(summary.scattered_dropped_by_reason);
        result["scattered_thinned_out"] = summary.scattered_thinned_out;
        return result;
      },
      nb::arg("input_path"), nb::arg("output_path"), nb::arg("transmission_path"),
      nb::arg("site_altitude_m"), nb::arg("profile_path"), nb::arg("envelope_paths"),
      nb::arg("wavelength_range"), nb::arg("seed"), nb::arg("mode"), nb::arg("store_p"),
      nb::arg("scattering_path"), nb::arg("scattering_points"), nb::arg("pathlen_m"),
      nb::arg("telescope_pathlen_m"), nb::arg("geometry_only"), nb::arg("curvature"),
      nb::arg("earth_radius_m"), nb::arg("telescopes"), nb::arg("threads"),
      nb::arg("scattering_primaries"), nb::arg("wavelength_interpolation"),
      nb::arg("dispersion_reference_nm"), nb::arg("ground_representation"),
      "Run the atmosphere stage on a photon stream file (geometry_only: no extinction, see "
      "docs/propagation/propagation.md). envelope_paths: (name, file) pairs of two-column text "
      "tables (nm, value). scattering_path: optional scattering table (single scattering of the "
      "primaries, off when None). pathlen_m: path inside each telescope [m]; telescope_pathlen_m: "
      "{telescope_id (1-based): metres} overrides. curvature: none, corsika or spherical (see "
      "Curvature); earth_radius_m: Earth radius [m], None for the CORSIKA CURVED radius. "
      "scattering_primaries: ids (1-based) whose photons are scattered (needs scattering; None: "
      "all blocks). telescopes: ids (1-based) to keep, None for all; other telescopes' blocks are "
      "dropped and not decoded. threads: 1 to 256 worker threads, the output does not depend on "
      "it. wavelength_interpolation: linear or floor (sim_telarray), see WavelengthInterpolation. "
      "dispersion_reference_nm: wavelength [nm] of the profile's n-1; None: no dispersion (needs "
      "a profile; rejected for refracted input). ground_representation: plane or sphere (rows on "
      "the telescope sphere, needs telescope_r_m; extinction of sphere rows is applied along the "
      "chord). The GIL is released while the stream is processed.");
}
