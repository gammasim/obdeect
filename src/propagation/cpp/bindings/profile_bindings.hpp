#pragma once

// Python API for the atmospheric profile, refraction tables and timing. Included by core.cpp
// after the nanobind headers. Arrays are one-dimensional float64 and are copied per call.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "dlufd/profile.hpp"
#include "dlufd/refraction.hpp"
#include "dlufd/timing.hpp"
#include "extinction_bindings.hpp" // to_numpy

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace dlufd_bindings {

using ConstArray = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

inline std::vector<double> copy_values(const ConstArray &values) {
  std::vector<double> copy(values.shape(0));
  for (std::size_t i = 0; i < copy.size(); ++i)
    copy[i] = values(i);
  return copy;
}

template <class Function>
nb::ndarray<nb::numpy, double, nb::ndim<1>> map_values(const ConstArray &heights,
                                                       Function function) {
  std::vector<double> out(heights.shape(0));
  for (std::size_t i = 0; i < out.size(); ++i)
    out[i] = function(heights(i));
  return to_numpy(std::move(out));
}

} // namespace dlufd_bindings

inline void bind_profile(nb::module_ &module) {
  using dlufd::AtmosphereProfile;
  using dlufd::RefractionTables;
  using dlufd_bindings::ConstArray;
  using dlufd_bindings::copy_values;
  using dlufd_bindings::map_values;
  using dlufd_bindings::to_numpy;

  nb::class_<AtmosphereProfile>(module, "AtmosphereProfile")
      .def(nb::new_([](const std::string &path) { return AtmosphereProfile::from_file(path); }),
           nb::arg("path"),
           "Read an atmprof table (altitude km, density, "
           "thickness, n-1).")
      .def_prop_ro("levels", &AtmosphereProfile::levels)
      .def_prop_ro("bottom_m", &AtmosphereProfile::bottom_m)
      .def_prop_ro("top_m", &AtmosphereProfile::top_m)
      .def(
          "density",
          [](const AtmosphereProfile &p, const ConstArray &h) {
            return map_values(h, [&p](double x) { return p.density(x); });
          },
          nb::arg("height_m"), "Density [g/cm^3] at heights [m a.s.l.] (spline table lookup).")
      .def(
          "thickness",
          [](const AtmosphereProfile &p, const ConstArray &h) {
            return map_values(h, [&p](double x) { return p.thickness(x); });
          },
          nb::arg("height_m"), "Thickness [g/cm^2] at heights [m a.s.l.] (spline table lookup).")
      .def(
          "refractive_index_minus_one",
          [](const AtmosphereProfile &p, const ConstArray &h) {
            return map_values(h, [&p](double x) { return p.refractive_index_minus_one(x); });
          },
          nb::arg("height_m"), "n-1 at heights [m a.s.l.] (spline table lookup).")
      .def(
          "group_index",
          [](const AtmosphereProfile &p, double observation_level_m) {
            return dlufd::ambient_group_index(p, observation_level_m);
          },
          nb::arg("observation_level_m"),
          "Refractive index at the observation level (obdeect ambient_group_index).");

  nb::class_<RefractionTables>(module, "RefractionTables")
      .def(nb::init<AtmosphereProfile, double>(), nb::arg("profile"),
           nb::arg("observation_level_m"),
           "Refraction tables for a profile and observation level [m a.s.l.].")
      .def_prop_ro("etadsn", &RefractionTables::etadsn)
      .def_prop_ro("observation_level_m",
                   [](const RefractionTables &t) { return t.observation_level_cm() / 100.; });

  module.def(
      "refract",
      [](const RefractionTables &tables, const ConstArray &h, const ConstArray &cx,
         const ConstArray &cy, const ConstArray &x, const ConstArray &y, const ConstArray &t) {
        const std::size_t n = h.shape(0);
        if (cx.shape(0) != n || cy.shape(0) != n || x.shape(0) != n || y.shape(0) != n ||
            t.shape(0) != n)
          throw std::invalid_argument("refract: all arrays must have the same length");
        std::vector<double> ux = copy_values(cx), uy = copy_values(cy), px = copy_values(x),
                            py = copy_values(y), pt = copy_values(t);
        std::vector<bool> bent(n);
        for (std::size_t i = 0; i < n; ++i) {
          const double height = h(i);
          bent[i] = dlufd::bend(tables, height, ux[i], uy[i], px[i], py[i], pt[i]);
        }
        std::vector<std::uint8_t> flags(bent.begin(), bent.end());
        return nb::make_tuple(to_numpy(std::move(ux)), to_numpy(std::move(uy)),
                              to_numpy(std::move(px)), to_numpy(std::move(py)),
                              to_numpy(std::move(pt)), to_numpy(std::move(flags)));
      },
      nb::arg("tables"), nb::arg("height_m"), nb::arg("cx"), nb::arg("cy"), nb::arg("x_m"),
      nb::arg("y_m"), nb::arg("t_ns"),
      "Bend photons emitted at height_m with direction (cx, cy) and position (x_m, y_m) at the "
      "observation level; returns (cx, cy, x_m, y_m, t_ns, bent) as numpy arrays, bent as uint8.");

  module.def("straight_path_air_delay_ns", &dlufd::straight_path_air_delay_ns, nb::arg("profile"),
             nb::arg("height_low_m"), nb::arg("height_high_m"), nb::arg("zenith_rad"),
             "Extra light travel time over a straight path between two heights [ns].");
  module.def("propagation_time_ns", &dlufd::propagation_time_ns, nb::arg("path_m"),
             nb::arg("group_index"), "Light travel time over a path in a medium [ns].");
}
