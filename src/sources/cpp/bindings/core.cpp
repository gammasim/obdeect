// Python.h must precede all standard headers (it sets feature-test macros such as
// _POSIX_C_SOURCE), so the nanobind headers come first.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include "dlamb/density.hpp"
#include "dlamb/source.hpp"
#include "dlamb/streams.hpp"
#include "dlamb/version.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

// Profiles cross the boundary boxed, so that the Python classes own a C++ value.
template <class Variant> struct Boxed {
  Variant value;
};
using SpectrumBox = Boxed<dlamb::Spectrum>;
using PulseBox = Boxed<dlamb::Pulse>;
using AngularBox = Boxed<dlamb::Angular>;

// Run a constructor that reports invalid input with std::invalid_argument as ValueError.
template <class F> auto checked(F &&f) -> decltype(f()) {
  try {
    return f();
  } catch (const std::invalid_argument &error) {
    throw nb::value_error(error.what());
  }
}

void need(const std::vector<double> &params, std::size_t count, const std::string &kind) {
  if (params.size() != count)
    throw nb::value_error(
        ("'" + kind + "' takes " + std::to_string(count) + " parameters").c_str());
}

dlamb::Vec3 vec3(const std::vector<double> &v, const char *name) {
  if (v.size() != 3)
    throw nb::value_error((std::string(name) + " must have three components").c_str());
  const dlamb::Vec3 result{v[0], v[1], v[2]};
  if (!dlamb::is_finite(result))
    throw nb::value_error((std::string(name) + " must be finite").c_str());
  return result;
}

// Angular profile. Parameters in degrees for the analytic kinds; `x` and `y` for "table"
// are u = 1 - cos(theta) and the intensity per solid angle.
// Half angle in degrees as u = 1 - cos(theta); the angle must be finite and in [0, 180].
double cone_u(const std::vector<double> &params, const std::string &kind) {
  const double degrees = params[0];
  if (!(degrees >= 0.0 && degrees <= 180.0))
    throw nb::value_error(("'" + kind + "' angle must be finite and in [0, 180] deg").c_str());
  return 1.0 - std::cos(degrees * dlamb::pi / 180.0);
}

AngularBox make_angular(const std::string &kind, const std::vector<double> &params,
                        const std::vector<double> &x, const std::vector<double> &y) {
  const double deg = dlamb::pi / 180.0;
  if (kind == "isotropic") {
    need(params, 0, kind);
    return {dlamb::Isotropic{}};
  }
  if (kind == "cone") {
    need(params, 1, kind);
    return {dlamb::FixedAngle{cone_u(params, kind)}};
  }
  if (kind == "filled_cone") {
    need(params, 1, kind);
    return {dlamb::FilledCone{cone_u(params, kind)}};
  }
  if (kind == "lambertian") {
    need(params, 0, kind);
    return {dlamb::Lambertian{}};
  }
  if (kind == "gauss") {
    need(params, 1, kind);
    return checked([&] {
      return AngularBox{dlamb::TabulatedAngular{dlamb::gaussian_angular_table(params[0] * deg)}};
    });
  }
  if (kind == "rayleigh") {
    need(params, 0, kind);
    return {dlamb::TabulatedAngular{dlamb::rayleigh_angular_table()}};
  }
  if (kind == "table") {
    need(params, 0, kind);
    return checked([&] { return AngularBox{dlamb::TabulatedAngular{{x, y}}}; });
  }
  throw nb::value_error(("unknown angular profile '" + kind + "'").c_str());
}

// Wavelength distribution [nm]. Parameters: fixed [nm]; uniform and inverse_square
// [low, high]; table: x = wavelength, y = density.
SpectrumBox make_spectrum(const std::string &kind, const std::vector<double> &params,
                          const std::vector<double> &x, const std::vector<double> &y) {
  if (kind == "fixed") {
    need(params, 1, kind);
    return {dlamb::Fixed{params[0]}};
  }
  if (kind == "uniform") {
    need(params, 2, kind);
    return {dlamb::UniformRange{params[0], params[1]}};
  }
  if (kind == "inverse_square") {
    need(params, 2, kind);
    return {dlamb::InverseSquare{params[0], params[1]}};
  }
  if (kind == "table") {
    need(params, 0, kind);
    return checked([&] { return SpectrumBox{dlamb::Tabulated{{x, y}}}; });
  }
  throw nb::value_error(("unknown spectrum '" + kind + "'").c_str());
}

// Emission time [ns]. Parameters: fixed [t]; tophat [full width, delay]; gauss [sigma,
// mean]; exponential [decay time, mean]; complex [flat width, gauss rms, decay time,
// delay]; table: x = time, y = rate.
PulseBox make_pulse(const std::string &kind, const std::vector<double> &params,
                    const std::vector<double> &x, const std::vector<double> &y) {
  if (kind == "fixed") {
    need(params, 1, kind);
    return {dlamb::Fixed{params[0]}};
  }
  if (kind == "tophat") {
    need(params, 2, kind);
    return {dlamb::UniformRange{params[1] - 0.5 * params[0], params[1] + 0.5 * params[0]}};
  }
  if (kind == "gauss") {
    need(params, 2, kind);
    return {dlamb::GaussianShape{params[1], params[0]}};
  }
  if (kind == "exponential") {
    need(params, 2, kind);
    return {dlamb::Exponential{params[1], params[0]}};
  }
  if (kind == "complex") {
    need(params, 4, kind);
    return {dlamb::ComplexPulse{params[0], params[1], params[2], params[3]}};
  }
  if (kind == "table") {
    need(params, 0, kind);
    return checked([&] { return PulseBox{dlamb::Tabulated{{x, y}}}; });
  }
  throw nb::value_error(("unknown pulse shape '" + kind + "'").c_str());
}

// NumPy array owning its buffer.
nb::object owned_float_array(std::vector<float> &&values) {
  auto *data = new std::vector<float>(std::move(values));
  nb::capsule owner(data, [](void *p) noexcept { delete static_cast<std::vector<float> *>(p); });
  return nb::ndarray<nb::numpy, float>(data->data(), {data->size()}, owner).cast();
}

nb::dict emit_chunk(const dlamb::Source &source, std::uint64_t seed, std::uint32_t event,
                    std::uint32_t telescope, std::uint64_t first_index, std::size_t rows) {
  gstrickd::PhotonChunk chunk;
  checked([&] {
    dlamb::emit(source, seed, event, telescope, first_index, rows, chunk);
    return 0;
  });
  nb::dict result;
  for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
    if ((dlamb::source_columns >> c) & 1U)
      result[gstrickd::column_name(static_cast<gstrickd::Column>(c))] =
          owned_float_array(std::move(chunk.columns[c]));
  }
  return result;
}

// Format constants of the gstrickd headers this build compiles against, and the columns
// dlamb emits. dlamb.source keeps Python copies of them; the tests compare the copies.
nb::dict gstrickd_contract() {
  nb::dict frames;
  frames["sky"] = static_cast<std::uint32_t>(gstrickd::Frame::sky);
  frames["telescope"] = static_cast<std::uint32_t>(gstrickd::Frame::telescope);
  frames["camera"] = static_cast<std::uint32_t>(gstrickd::Frame::camera);
  std::vector<std::string> emitted;
  for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
    if ((dlamb::source_columns >> c) & 1U)
      emitted.emplace_back(gstrickd::column_name(static_cast<gstrickd::Column>(c)));
  }
  nb::dict result;
  result["format_major"] = gstrickd::format_major;
  result["format_minor"] = gstrickd::format_minor;
  result["max_chunk_rows"] = gstrickd::max_chunk_rows;
  result["dlamb_format_major"] = dlamb::photon_stream_format_major;
  result["frames"] = frames;
  result["emitted_columns"] = emitted;
  return result;
}

} // namespace

NB_MODULE(_core, module) {
  module.doc() = "dlamb native core";
  module.def("core_version", [] {
    return std::make_tuple(dlamb::version_major, dlamb::version_minor, dlamb::version_patch);
  });
  module.def("gstrickd_contract", &gstrickd_contract,
             "Format constants of the gstrickd headers and the columns dlamb emits.");

  nb::class_<SpectrumBox>(module, "Spectrum");
  nb::class_<PulseBox>(module, "Pulse");
  nb::class_<AngularBox>(module, "Angular");
  module.def("make_spectrum", &make_spectrum, "kind"_a, "params"_a, "x"_a, "y"_a);
  module.def("make_pulse", &make_pulse, "kind"_a, "params"_a, "x"_a, "y"_a);
  module.def("make_angular", &make_angular, "kind"_a, "params"_a, "x"_a, "y"_a);

  nb::class_<dlamb::Emitter>(module, "Emitter")
      .def(nb::new_([](const std::vector<double> &position, const std::vector<double> &pointing,
                       const AngularBox &angular) {
             return dlamb::Emitter{vec3(position, "position"), vec3(pointing, "pointing"),
                                   angular.value};
           }),
           "position"_a, "pointing"_a, "angular"_a);

  nb::class_<dlamb::Star>(module, "Star")
      .def(nb::new_([](const std::vector<double> &towards, double distance_m, double disk_radius_m,
                       double disk_z_m, double start_z_m) {
             return dlamb::Star{vec3(towards, "towards_source"), distance_m, disk_radius_m,
                                disk_z_m, start_z_m};
           }),
           "towards_source"_a, "distance_m"_a, "disk_radius_m"_a, "disk_z_m"_a, "start_z_m"_a);

  nb::class_<dlamb::Source>(module, "Source")
      .def(nb::new_([](int frame, const nb::object &geometry, const SpectrumBox &spectrum,
                       const PulseBox &pulse, double photons) {
             dlamb::Source source;
             if (frame < 0 || frame > 3)
               throw nb::value_error("frame must be 0 (sky), 2 (T) or 3 (C)");
             source.frame = static_cast<gstrickd::Frame>(frame);
             if (nb::isinstance<dlamb::Emitter>(geometry))
               source.geometry = nb::cast<dlamb::Emitter>(geometry);
             else if (nb::isinstance<dlamb::Star>(geometry))
               source.geometry = nb::cast<dlamb::Star>(geometry);
             else
               throw nb::type_error("geometry must be an Emitter or a Star");
             source.spectrum = spectrum.value;
             source.pulse = pulse.value;
             source.photons = photons;
             checked([&] {
               dlamb::validate(source);
               return 0;
             });
             return source;
           }),
           "frame"_a, "geometry"_a, "spectrum"_a, "pulse"_a, "photons"_a)
      .def(
          "set_disc_target",
          [](dlamb::Source &self, const std::vector<double> &centre,
             const std::vector<double> &normal, double radius_m) {
            self.target =
                dlamb::DiscTarget{vec3(centre, "centre"), vec3(normal, "normal"), radius_m};
          },
          "centre"_a, "normal"_a, "radius_m"_a)
      .def(
          "set_sphere_target",
          [](dlamb::Source &self, const std::vector<double> &centre, double radius_m) {
            self.target = dlamb::SphereTarget{vec3(centre, "centre"), radius_m};
          },
          "centre"_a, "radius_m"_a)
      .def("validate",
           [](const dlamb::Source &self) {
             checked([&] {
               dlamb::validate(self);
               return 0;
             });
           })
      .def("expected_count",
           [](const dlamb::Source &self) {
             return checked([&] { return dlamb::expected_count(self); });
           })
      .def(
          "photon_count",
          [](const dlamb::Source &self, std::uint64_t seed, std::uint32_t event,
             std::uint32_t telescope) {
            return checked([&] { return dlamb::photon_count(self, seed, event, telescope); });
          },
          "seed"_a, "event"_a, "telescope"_a)
      .def("emit", &emit_chunk, "seed"_a, "event"_a, "telescope"_a, "first_index"_a, "rows"_a);
}
