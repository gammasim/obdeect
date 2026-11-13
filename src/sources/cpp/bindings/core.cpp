// nanobind extension dlamb._core: the only C++ file that depends on Python.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "obdeect/sources/density.hpp"
#include "obdeect/sources/source.hpp"
#include "obdeect/sources/streams.hpp"
#include "obdeect/sources/version.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;
using namespace obdeect::sources;

namespace {

Vec3 to_vec3(const std::vector<double> &v, const char *what) {
  if (v.size() != 3)
    throw std::invalid_argument(std::string(what) + " needs three components");
  const Vec3 out{v[0], v[1], v[2]};
  if (!is_finite(out))
    throw std::invalid_argument(std::string(what) + " must be finite");
  return out;
}

void need(const std::vector<double> &params, std::size_t n, const char *kind) {
  if (params.size() != n)
    throw std::invalid_argument(std::string(kind) + ": wrong number of parameters");
}

LinearDensity table(const std::vector<double> &x, const std::vector<double> &y) {
  return LinearDensity(x, y);
}

Spectrum make_spectrum(const std::string &kind, const std::vector<double> &p,
                       const std::vector<double> &x, const std::vector<double> &y) {
  if (kind == "fixed") {
    need(p, 1, "fixed spectrum");
    return Fixed{p[0]};
  }
  if (kind == "uniform") {
    need(p, 2, "uniform spectrum");
    return UniformRange{p[0], p[1]};
  }
  if (kind == "inverse_square") {
    need(p, 2, "inverse_square spectrum");
    return InverseSquare{p[0], p[1]};
  }
  if (kind == "table") {
    need(p, 0, "table spectrum");
    return Tabulated{table(x, y)};
  }
  throw std::invalid_argument("unknown spectrum kind '" + kind + "'");
}

Pulse make_pulse(const std::string &kind, const std::vector<double> &p,
                 const std::vector<double> &x, const std::vector<double> &y) {
  if (kind == "fixed") {
    need(p, 1, "fixed pulse");
    return Fixed{p[0]};
  }
  if (kind == "tophat") {
    need(p, 2, "tophat pulse");
    return UniformRange{p[1] - 0.5 * p[0], p[1] + 0.5 * p[0]};
  }
  if (kind == "gauss") {
    need(p, 2, "gauss pulse");
    return GaussianShape{p[1], p[0]};
  }
  if (kind == "exponential") {
    need(p, 2, "exponential pulse");
    return Exponential{p[1], p[0]};
  }
  if (kind == "complex") {
    need(p, 4, "complex pulse");
    return ComplexPulse{p[0], p[1], p[2], p[3]};
  }
  if (kind == "table") {
    need(p, 0, "table pulse");
    return Tabulated{table(x, y)};
  }
  throw std::invalid_argument("unknown pulse kind '" + kind + "'");
}

double cone_u0(double angle_deg) {
  if (!(std::isfinite(angle_deg) && angle_deg >= 0.0 && angle_deg <= 180.0))
    throw std::invalid_argument("angle_deg must be finite and in [0, 180]");
  return 1.0 - std::cos(angle_deg * pi / 180.0);
}

Angular make_angular(const std::string &kind, const std::vector<double> &p,
                     const std::vector<double> &x, const std::vector<double> &y) {
  if (kind == "isotropic" || kind == "lambertian" || kind == "rayleigh" || kind == "table")
    need(p, 0, kind.c_str());
  if (kind == "isotropic")
    return Isotropic{};
  if (kind == "lambertian")
    return Lambertian{};
  if (kind == "rayleigh")
    return TabulatedAngular{rayleigh_angular_table()};
  if (kind == "cone") {
    need(p, 1, "cone");
    return FixedAngle{cone_u0(p[0])};
  }
  if (kind == "filled_cone") {
    need(p, 1, "filled_cone");
    return FilledCone{cone_u0(p[0])};
  }
  if (kind == "gauss") {
    need(p, 1, "gauss");
    return TabulatedAngular{gaussian_angular_table(p[0] * pi / 180.0)};
  }
  if (kind == "table")
    return TabulatedAngular{table(x, y)};
  throw std::invalid_argument("unknown angular kind '" + kind + "'");
}

// A float32 numpy array that owns a copy of the values.
nb::object to_array(const std::vector<float> &values) {
  auto *data = new float[values.size()];
  std::copy(values.begin(), values.end(), data);
  nb::capsule owner(data, [](void *p) noexcept { delete[] static_cast<float *>(p); });
  return nb::cast(nb::ndarray<nb::numpy, float>(data, {values.size()}, owner));
}

} // namespace

NB_MODULE(_core, m) {
  m.doc() = "dlamb native core";
  m.def("core_version", [] { return nb::make_tuple(version_major, version_minor, version_patch); });
  m.def("stream_contract", [] {
    nb::dict d;
    d["format_major"] = obdeect::stream::format_major;
    d["format_minor"] = obdeect::stream::format_minor;
    d["max_chunk_rows"] = obdeect::stream::max_chunk_rows;
    d["sources_format_major"] = photon_stream_format_major;
    nb::dict frames;
    frames["sky"] = static_cast<int>(obdeect::stream::Frame::sky);
    frames["telescope"] = static_cast<int>(obdeect::stream::Frame::telescope);
    frames["camera"] = static_cast<int>(obdeect::stream::Frame::camera);
    d["frames"] = frames;
    nb::list columns;
    for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c)
      if ((source_columns >> c) & 1U)
        columns.append(obdeect::stream::column_name(static_cast<obdeect::stream::Column>(c)));
    d["emitted_columns"] = columns;
    return d;
  });

  nb::class_<Spectrum>(m, "Spectrum");
  nb::class_<Pulse>(m, "Pulse");
  nb::class_<Angular>(m, "Angular");
  m.def("make_spectrum", &make_spectrum, nb::arg("kind"), nb::arg("params"), nb::arg("x"),
        nb::arg("y"));
  m.def("make_pulse", &make_pulse, nb::arg("kind"), nb::arg("params"), nb::arg("x"), nb::arg("y"));
  m.def("make_angular", &make_angular, nb::arg("kind"), nb::arg("params"), nb::arg("x"),
        nb::arg("y"));

  nb::class_<Emitter>(m, "Emitter")
      .def("__init__", [](Emitter *self, const std::vector<double> &position,
                          const std::vector<double> &pointing, const Angular &angular) {
        new (self) Emitter{to_vec3(position, "position"), to_vec3(pointing, "pointing"), angular};
      });
  nb::class_<Star>(m, "Star").def("__init__", [](Star *self, const std::vector<double> &towards,
                                                 double distance_m, double disk_radius_m,
                                                 double disk_z_m, double start_z_m) {
    new (self)
        Star{to_vec3(towards, "towards_source"), distance_m, disk_radius_m, disk_z_m, start_z_m};
  });

  nb::class_<Source>(m, "Source")
      .def("__init__",
           [](Source *self, int frame, nb::object geometry, const Spectrum &spectrum,
              const Pulse &pulse, double photons) {
             if (frame != 0 && frame != 2 && frame != 3)
               throw std::invalid_argument("frame must be 0 (sky), 2 (telescope) or 3 (camera)");
             new (self) Source();
             self->frame = static_cast<obdeect::stream::Frame>(frame);
             if (nb::isinstance<Emitter>(geometry))
               self->geometry = nb::cast<Emitter>(geometry);
             else if (nb::isinstance<Star>(geometry))
               self->geometry = nb::cast<Star>(geometry);
             else
               throw std::invalid_argument("geometry must be an Emitter or a Star");
             self->spectrum = spectrum;
             self->pulse = pulse;
             self->photons = photons;
             validate(*self);
           })
      .def("set_disc_target",
           [](Source &s, const std::vector<double> &centre, const std::vector<double> &normal,
              double radius) {
             s.target = DiscTarget{to_vec3(centre, "centre"), to_vec3(normal, "normal"), radius};
           })
      .def("set_sphere_target",
           [](Source &s, const std::vector<double> &centre, double radius) {
             s.target = SphereTarget{to_vec3(centre, "centre"), radius};
           })
      .def("validate", [](const Source &s) { validate(s); })
      .def("expected_count", [](const Source &s) { return expected_count(s); })
      .def("photon_count",
           [](const Source &s, std::uint64_t seed, std::uint32_t event, std::uint32_t telescope) {
             return photon_count(s, seed, event, telescope);
           })
      .def("emit", [](const Source &s, std::uint64_t seed, std::uint32_t event,
                      std::uint32_t telescope, std::uint64_t first_index, std::size_t rows) {
        obdeect::stream::PhotonChunk chunk;
        emit(s, seed, event, telescope, first_index, rows, chunk);
        nb::dict out;
        for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c)
          if ((source_columns >> c) & 1U) {
            const auto column = static_cast<obdeect::stream::Column>(c);
            out[obdeect::stream::column_name(column)] = to_array(chunk[column]);
          }
        return out;
      });
}
