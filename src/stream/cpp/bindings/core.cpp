// Python.h must precede all standard headers (it sets feature-test macros such as
// _POSIX_C_SOURCE), so the nanobind headers come first.
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

#include "gstrickd/version.hpp"

#include <tuple>

namespace nb = nanobind;

NB_MODULE(_core, module) {
  module.doc() = "gstrickd native core";
  module.def("core_version", [] {
    return std::make_tuple(gstrickd::version_major, gstrickd::version_minor,
                           gstrickd::version_patch);
  });
}
