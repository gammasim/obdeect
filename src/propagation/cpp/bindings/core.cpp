// Python.h must precede all standard headers (it sets feature-test macros such as
// _POSIX_C_SOURCE), so the nanobind headers come first.
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

#include "obdeect/propagation/version.hpp"
#include "extinction_bindings.hpp"
#include "profile_bindings.hpp"
#include "propagator_bindings.hpp"
#include "scattering_bindings.hpp"

#include <tuple>

namespace nb = nanobind;

NB_MODULE(_core, module) {
  module.doc() = "dlufd native core";
  module.def("core_version", [] {
    return std::make_tuple(obdeect::propagation::version_major, obdeect::propagation::version_minor, obdeect::propagation::version_patch);
  });
  bind_extinction(module);
  bind_profile(module);
  bind_propagator(module);
  bind_scattering(module);
}
