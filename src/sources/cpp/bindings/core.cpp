// Python.h must precede all standard headers (it sets feature-test macros such as
// _POSIX_C_SOURCE), so the nanobind headers come first.
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

#include "dlamb/version.hpp"

#include <tuple>

namespace nb = nanobind;

NB_MODULE(_core, module) {
  module.doc() = "dlamb native core";
  module.def("core_version", [] {
    return std::make_tuple(dlamb::version_major, dlamb::version_minor, dlamb::version_patch);
  });
}
