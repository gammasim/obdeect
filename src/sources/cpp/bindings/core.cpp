// nanobind extension obdeect.sources._core: the only C++ file that depends on Python.
//
// The Python side passes normalized dictionaries (obdeect.sources.config.normalize): a run is
// {"run": {...}, "telescopes": [...], "events": [...]}, an emission is one dictionary with a
// "kind" (point, line, plane_wave, sky_background). Every dictionary is converted and validated
// again here, with the key path in the error message.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "obdeect/sources/density.hpp"
#include "obdeect/sources/emit.hpp"
#include "obdeect/sources/run.hpp"
#include "obdeect/sources/source.hpp"
#include "obdeect/sources/streams.hpp"
#include "obdeect/sources/version.hpp"

#include "obdeect/stream/format.hpp"
#include "obdeect/stream/photons.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;
using namespace obdeect::sources;

namespace {

[[noreturn]] void bad(const std::string &where, const std::string &what) {
  throw std::invalid_argument(where + ": " + what);
}

// A dictionary of a normalized description, with its key path for error messages.
class Section {
public:
  Section(nb::handle h, std::string path) : path_(std::move(path)) {
    if (!nb::isinstance<nb::dict>(h))
      bad(path_, "expected a dictionary");
    dict_ = nb::borrow<nb::dict>(h);
  }

  bool has(const char *key) const { return dict_.contains(nb::str(key)); }

  std::string where(const char *key) const {
    return path_.empty() ? std::string(key) : path_ + "." + key;
  }

  nb::object get(const char *key) const {
    if (!has(key))
      bad(where(key), "missing");
    nb::object value = dict_[nb::str(key)];
    return value;
  }

  bool has_value(const char *key) const { return has(key) && !get(key).is_none(); }

  double real(const char *key) const {
    const nb::object v = get(key);
    if (nb::isinstance<nb::bool_>(v) ||
        !(nb::isinstance<nb::int_>(v) || nb::isinstance<nb::float_>(v)))
      bad(where(key), "expected a number");
    double out{};
    if (!nb::try_cast(v, out) || !std::isfinite(out))
      bad(where(key), "expected a finite number");
    return out;
  }

  std::int64_t integer(const char *key) const {
    const nb::object v = get(key);
    std::int64_t out{};
    if (nb::isinstance<nb::bool_>(v) || !nb::isinstance<nb::int_>(v) || !nb::try_cast(v, out))
      bad(where(key), "expected an integer");
    return out;
  }

  std::uint64_t unsigned_integer(const char *key) const {
    const nb::object v = get(key);
    std::uint64_t out{};
    if (nb::isinstance<nb::bool_>(v) || !nb::isinstance<nb::int_>(v) || !nb::try_cast(v, out))
      bad(where(key), "expected an integer from 0 to 2^64 - 1");
    return out;
  }

  bool flag(const char *key) const {
    const nb::object v = get(key);
    if (!nb::isinstance<nb::bool_>(v))
      bad(where(key), "expected true or false");
    return nb::cast<bool>(v);
  }

  std::string text(const char *key) const {
    const nb::object v = get(key);
    std::string out;
    if (!nb::isinstance<nb::str>(v) || !nb::try_cast(v, out))
      bad(where(key), "expected a string");
    return out;
  }

  std::vector<double> numbers(const char *key) const {
    const nb::object v = get(key);
    std::vector<double> out;
    if (nb::isinstance<nb::str>(v) || !nb::try_cast(v, out))
      bad(where(key), "expected a list of numbers");
    for (const double x : out)
      if (!std::isfinite(x))
        bad(where(key), "values must be finite");
    return out;
  }

  Vec3 vec3(const char *key) const {
    const auto v = numbers(key);
    if (v.size() != 3)
      bad(where(key), "needs three components");
    return Vec3{v[0], v[1], v[2]};
  }

  // A unit vector (the direction is normalised; a zero vector is an error).
  Vec3 unit(const char *key) const {
    const Vec3 v = vec3(key);
    const double n = norm(v);
    if (!(n > 0.0))
      bad(where(key), "needs a non-zero vector");
    return v * (1.0 / n);
  }

  Section section(const char *key) const { return Section(get(key), where(key)); }

  // Items of a list; `f(Section)` is called for each with its path "key[i]".
  template <class F> void each(const char *key, F &&f) const {
    const nb::object v = get(key);
    if (!nb::isinstance<nb::list>(v))
      bad(where(key), "expected a list");
    std::size_t i = 0;
    for (nb::handle item : nb::borrow<nb::list>(v)) {
      f(Section(item, where(key) + "[" + std::to_string(i) + "]"));
      ++i;
    }
  }

private:
  nb::dict dict_;
  std::string path_;
};

obdeect::stream::Value value_of(nb::handle v, const std::string &where) {
  if (nb::isinstance<nb::bool_>(v))
    bad(where, "booleans are not metadata values");
  if (nb::isinstance<nb::int_>(v)) {
    std::int64_t x{};
    if (!nb::try_cast(v, x))
      bad(where, "integer out of range");
    return x;
  }
  if (nb::isinstance<nb::float_>(v)) {
    double x{};
    nb::try_cast(v, x);
    return x;
  }
  if (nb::isinstance<nb::str>(v)) {
    std::string x;
    nb::try_cast(v, x);
    return x;
  }
  std::vector<double> x;
  if (nb::try_cast(v, x))
    return x;
  bad(where, "metadata values must be int, float, str or a list of floats");
}

obdeect::stream::Metadata metadata_of(const Section &s, const char *key) {
  obdeect::stream::Metadata out;
  if (!s.has_value(key))
    return out;
  const nb::object v = s.get(key);
  if (!nb::isinstance<nb::dict>(v))
    bad(s.where(key), "expected a dictionary");
  for (auto [k, value] : nb::borrow<nb::dict>(v)) {
    if (!nb::isinstance<nb::str>(k))
      bad(s.where(key), "metadata keys must be strings");
    std::string name;
    nb::try_cast(k, name);
    out.emplace(name, value_of(value, s.where(key) + "." + name));
  }
  return out;
}

obdeect::stream::Frame frame_of(const std::string &name, const std::string &where) {
  if (name == "sky")
    return obdeect::stream::Frame::sky;
  if (name == "telescope")
    return obdeect::stream::Frame::telescope;
  if (name == "camera")
    return obdeect::stream::Frame::camera;
  bad(where, "frame must be sky, telescope or camera, not '" + name + "'");
}

double cone_u0(double angle_deg, const std::string &where) {
  if (!(std::isfinite(angle_deg) && angle_deg >= 0.0 && angle_deg <= 180.0))
    bad(where, "angle_deg must be finite and in [0, 180]");
  return 1.0 - std::cos(angle_deg * pi / 180.0);
}

Angular angular_of(const Section &s) {
  const std::string kind = s.text("kind");
  if (kind == "isotropic")
    return Isotropic{};
  if (kind == "lambertian")
    return Lambertian{};
  if (kind == "rayleigh")
    return TabulatedAngular{rayleigh_angular_table()};
  if (kind == "cone")
    return FixedAngle{cone_u0(s.real("angle_deg"), s.where("angle_deg"))};
  if (kind == "filled_cone")
    return FilledCone{cone_u0(s.real("angle_deg"), s.where("angle_deg"))};
  if (kind == "gauss") {
    const double sigma = s.real("sigma_deg");
    if (!(sigma >= 0.0 && sigma <= 180.0))
      bad(s.where("sigma_deg"), "sigma_deg must be in [0, 180]");
    return TabulatedAngular{gaussian_angular_table(sigma * pi / 180.0)};
  }
  if (kind == "table") {
    // The density per solid angle is tabulated in u = 1 - cos(theta).
    const auto degrees = s.numbers("angle_deg");
    std::vector<double> u;
    for (const double d : degrees) {
      if (!(d >= 0.0 && d <= 180.0))
        bad(s.where("angle_deg"), "angle_deg values must be in [0, 180]");
      u.push_back(u_of_angle(d * pi / 180.0));
    }
    return TabulatedAngular{LinearDensity(u, s.numbers("value"))};
  }
  bad(s.where("kind"), "unknown angular kind '" + kind + "'");
}

Spectrum spectrum_of(const Section &s) {
  const std::string kind = s.text("kind");
  if (kind == "fixed")
    return Fixed{s.real("wavelength_nm")};
  if (kind == "uniform")
    return UniformRange{s.real("low_nm"), s.real("high_nm")};
  if (kind == "inverse_square")
    return InverseSquare{s.real("low_nm"), s.real("high_nm")};
  if (kind == "table")
    return Tabulated{LinearDensity(s.numbers("wavelength_nm"), s.numbers("density"))};
  bad(s.where("kind"), "unknown spectrum kind '" + kind + "'");
}

Pulse pulse_of(const Section &s) {
  const std::string kind = s.text("kind");
  if (kind == "fixed")
    return Fixed{s.real("time_ns")};
  if (kind == "tophat") {
    const double width = s.real("width_ns");
    const double delay = s.real("delay_ns");
    return UniformRange{delay - 0.5 * width, delay + 0.5 * width};
  }
  if (kind == "gauss")
    return GaussianShape{s.real("mean_ns"), s.real("sigma_ns")};
  if (kind == "exponential")
    return Exponential{s.real("mean_ns"), s.real("decay_ns")};
  if (kind == "complex")
    return ComplexPulse{s.real("flat_ns"), s.real("gauss_rms_ns"), s.real("decay_ns"),
                        s.real("delay_ns")};
  if (kind == "table")
    return Tabulated{LinearDensity(s.numbers("time_ns"), s.numbers("rate"))};
  bad(s.where("kind"), "unknown pulse kind '" + kind + "'");
}

Target target_of(const Section &s) {
  const std::string kind = s.text("kind");
  if (kind == "none")
    return std::monostate{};
  if (kind == "auto")
    return AutoTarget{};
  if (kind == "disc") {
    DiscTarget d;
    d.centre = s.vec3("centre_m");
    d.normal = s.unit("normal");
    d.radius_m = s.real("radius_m");
    return d;
  }
  if (kind == "sphere") {
    SphereTarget t;
    t.centre = s.vec3("centre_m");
    t.radius_m = s.real("radius_m");
    return t;
  }
  bad(s.where("kind"), "unknown target kind '" + kind + "'");
}

Geometry geometry_of(const Section &s, const std::string &kind) {
  if (kind == "point") {
    Point p;
    p.position = s.vec3("position_m");
    p.axis = s.unit("axis");
    p.angular = angular_of(s.section("angular"));
    return p;
  }
  if (kind == "line") {
    Line l;
    l.start = s.vec3("start_m");
    l.direction = s.unit("direction");
    l.length_m = s.real("length_m");
    l.beta = s.real("beta");
    l.angular = angular_of(s.section("angular"));
    if (s.has_value("profile")) {
      const Section profile = s.section("profile");
      l.profile = LinearDensity(profile.numbers("distance_m"), profile.numbers("density"));
    }
    return l;
  }
  if (kind == "plane_wave") {
    PlaneWave w;
    w.towards_source = s.unit("towards_source");
    if (s.has_value("distance_m"))
      w.distance_m = s.real("distance_m");
    w.disc_radius_m = s.real("disc_radius_m");
    w.disc_z_m = s.real("disc_z_m");
    w.start_z_m = s.real("start_z_m");
    return w;
  }
  if (kind == "sky_background") {
    SkyBackground b;
    b.centre = s.vec3("centre_m");
    b.radius_m = s.real("radius_m");
    b.axis = s.unit("axis");
    b.viewcone_deg = s.real("viewcone_deg");
    return b;
  }
  bad(s.where("kind"), "unknown emission kind '" + kind + "'");
}

// An emission of a run: photons are the physical count (plan() divides by the bunch size).
Emission emission_of(const Section &s) {
  Emission e;
  e.source.geometry = geometry_of(s, s.text("kind"));
  e.source.spectrum = spectrum_of(s.section("spectrum"));
  e.source.pulse = pulse_of(s.section("pulse"));
  e.source.target = target_of(s.section("target"));
  e.source.photons = s.real("photons");
  e.bunch_size = s.real("bunch_size");
  e.time_offset_ns = s.real("time_offset_ns");
  e.to_plane_z0 = s.flag("to_plane_z0");
  return e;
}

// The emission of a single-emission call (no run): photons per bunch, as plan() does it.
Emission single_emission(nb::handle h) {
  const Section s(h, "emission");
  Emission e = emission_of(s);
  if (std::holds_alternative<AutoTarget>(e.source.target))
    bad("emission.target", "'auto' needs the telescopes of a run (use write or check)");
  if (!(std::isfinite(e.bunch_size) && e.bunch_size >= 1.0))
    bad("emission.bunch_size", "must be finite and >= 1");
  if (!within_limits(e.time_offset_ns))
    bad("emission.time_offset_ns", "out of range");
  e.source.photons /= e.bunch_size;
  validate(e.source);
  return e;
}

Run run_of(nb::handle h) {
  const Section top(h, "");
  Run run;
  const Section r = top.section("run");
  run.number = r.integer("number");
  run.seed = r.unsigned_integer("seed");
  run.observation_level_m = r.real("observation_level_m");
  run.atmosphere_id = r.integer("atmosphere_id");
  run.frame = frame_of(r.text("frame"), r.where("frame"));
  run.metadata = metadata_of(r, "metadata");
  top.each("telescopes", [&](const Section &t) {
    Telescope tel;
    const auto id = t.integer("id");
    if (id < 0 || id > std::numeric_limits<std::uint32_t>::max())
      bad(t.where("id"), "telescope id out of range");
    tel.id = static_cast<std::uint32_t>(id);
    tel.position = t.vec3("position_m");
    tel.radius_m = t.real("radius_m");
    run.telescopes.push_back(tel);
  });
  top.each("events", [&](const Section &ev) {
    Event event;
    event.repeat = ev.unsigned_integer("repeat");
    event.metadata = metadata_of(ev, "metadata");
    ev.each("emissions", [&](const Section &em) { event.emissions.push_back(emission_of(em)); });
    run.events.push_back(std::move(event));
  });
  return run;
}

// A float32 numpy array that owns a copy of the values.
nb::object to_array(const std::vector<float> &values) {
  auto *data = new float[values.size()];
  std::copy(values.begin(), values.end(), data);
  nb::capsule owner(data, [](void *p) noexcept { delete[] static_cast<float *>(p); });
  return nb::cast(nb::ndarray<nb::numpy, float>(data, {values.size()}, owner));
}

nb::dict write_run_py(nb::handle run_description, const std::string &path) {
  const Run run = run_of(run_description);
  RunSummary summary;
  {
    nb::gil_scoped_release release;
    if (path == "-") {
      summary = write_run(std::cout, run);
      std::cout.flush();
    } else {
      std::ofstream out(path, std::ios::binary);
      if (!out)
        throw std::runtime_error("cannot open " + path + " for writing");
      summary = write_run(out, run);
      out.flush();
      if (!out)
        throw std::runtime_error("writing " + path + " failed");
    }
  }
  nb::dict out;
  out["events"] = summary.events;
  out["rows"] = summary.rows;
  out["weight"] = summary.weight;
  out["expected"] = summary.expected;
  out["not_on_plane"] = summary.not_on_plane;
  return out;
}

nb::dict check_run_py(nb::handle run_description) {
  const Run run = run_of(run_description);
  const auto planned = plan(run);
  nb::list emissions;
  std::uint64_t events = 0;
  double candidates = 0.0;
  for (std::size_t i = 0; i < planned.size(); ++i) {
    const double repeat = static_cast<double>(run.events[i].repeat);
    events += run.events[i].repeat;
    for (std::size_t j = 0; j < planned[i].size(); ++j) {
      const auto &p = planned[i][j];
      const double mean = expected_count(p.emission.source);
      candidates += mean * repeat * static_cast<double>(run.telescopes.size());
      nb::dict d;
      d["event"] = i;
      d["emission"] = j;
      d["repeat"] = run.events[i].repeat;
      d["seed"] = p.seed;
      d["candidates"] = mean;
      d["expected_weight"] = expected_weight(p.emission.source) * p.emission.bunch_size;
      emissions.append(d);
    }
  }
  nb::dict out;
  out["events"] = events;
  out["telescopes"] = run.telescopes.size();
  out["emissions"] = emissions;
  out["expected_candidates"] = candidates;
  out["expected_weight"] = expected_weight(run);
  return out;
}

nb::dict emit_chunk_py(nb::handle emission, std::uint64_t seed, std::uint32_t event,
                       std::uint32_t telescope, std::uint64_t first_index, std::size_t rows) {
  const Emission e = single_emission(emission);
  obdeect::stream::PhotonChunk chunk;
  emit(e.source, seed, event, telescope, first_index, rows, chunk);
  // The run's per-emission rows (bunch weight, time offset, move to z = 0); no zem column. Rows
  // that do not reach the plane are dropped (their count is not returned here).
  Run camera_run;
  camera_run.frame = obdeect::stream::Frame::camera;
  std::uint64_t not_on_plane = 0;
  detail::finish_rows(e, camera_run, chunk, not_on_plane);
  nb::dict out;
  for (std::uint32_t c = 0; c < obdeect::stream::n_columns; ++c)
    if ((source_columns >> c) & 1U) {
      const auto column = static_cast<obdeect::stream::Column>(c);
      out[obdeect::stream::column_name(column)] = to_array(chunk[column]);
    }
  return out;
}

double expected_count_py(nb::handle emission) {
  return expected_count(single_emission(emission).source);
}

double expected_weight_py(nb::handle emission) {
  const Emission e = single_emission(emission);
  return expected_weight(e.source) * e.bunch_size;
}

std::uint64_t photon_count_py(nb::handle emission, std::uint64_t seed, std::uint32_t event,
                              std::uint32_t telescope) {
  return photon_count(single_emission(emission).source, seed, event, telescope);
}

} // namespace

NB_MODULE(_core, m) {
  m.doc() = "obdeect sources native core";
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
  m.def("write_run", &write_run_py, nb::arg("run"), nb::arg("path"),
        "Write a normalized run to a photon-stream file ('-' = stdout); returns the summary.");
  m.def("check_run", &check_run_py, nb::arg("run"),
        "Candidates and expected weights of a normalized run, without generating photons.");
  m.def("emit_chunk", &emit_chunk_py, nb::arg("emission"), nb::arg("seed"), nb::arg("event"),
        nb::arg("telescope"), nb::arg("first_index"), nb::arg("rows"));
  m.def("expected_count", &expected_count_py, nb::arg("emission"));
  m.def("expected_weight", &expected_weight_py, nb::arg("emission"));
  m.def("photon_count", &photon_count_py, nb::arg("emission"), nb::arg("seed"), nb::arg("event"),
        nb::arg("telescope"));
}
