#pragma once

// Trace a photon stream through optical models. One implementation for the native program
// `obdeect-ray-tracing-trace-stream` and the Python API (`_core.trace_stream`): ground-frame (G)
// or telescope-frame (T) photons in, camera-frame (C) photons out. Records pass through,
// telescopes not selected are dropped and the efficiency status is updated
// (docs/ray-tracing/photon-stream.md).

#include "obdeect/ray_tracing/bunches.hpp"
#include "obdeect/ray_tracing/json.hpp"
#include "obdeect/ray_tracing/model_hash.hpp"
#include "obdeect/ray_tracing/model_json.hpp"
#include "obdeect/stream/photons.hpp"
#include "obdeect/stream/random.hpp"
#include "obdeect/stream/stream.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <istream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace obdeect::ray_tracing {

namespace stream_trace_detail {

namespace st = obdeect::stream;
using json::Value;

inline constexpr double kPi = 3.141592653589793;
// Margin above the highest incoming obscurer for the input plane [m] (as simtools).
inline constexpr double kInputPlaneMarginM = 1.0;
// Radial samples of an asphere mirror's aperture for its support.
inline constexpr std::size_t kAsphereSamples = 2049;
inline constexpr const char *kConstantFactor = "telescope_transmission_constant";

using Vec = std::array<double, 3>;

[[nodiscard]] inline double dot(const Vec &a, const Vec &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// numpy.interp-like linear interpolation of (nodes, values) at x, nodes increasing, x inside
// [nodes.front(), nodes.back()] (checked by the caller).
[[nodiscard]] inline double interpolate(double x, const std::vector<double> &nodes,
                                        const std::vector<double> &values) {
  const std::size_t n = nodes.size();
  if (n == 1)
    return values[0];
  const auto above = std::upper_bound(nodes.begin(), nodes.end(), x);
  const std::size_t j = static_cast<std::size_t>(above - nodes.begin()) - 1;
  if (j == n - 1 || nodes[j] == x)
    return values[j];
  const double slope = (values[j + 1] - values[j]) / (nodes[j + 1] - nodes[j]);
  return slope * (x - nodes[j]) + values[j];
}

[[nodiscard]] inline std::string format_g(double value) {
  char buffer[64];
  std::snprintf(buffer, sizeof buffer, "%g", value);
  return buffer;
}

[[nodiscard]] inline std::string format_3f(double value) {
  char buffer[64];
  std::snprintf(buffer, sizeof buffer, "%.3f", value);
  return buffer;
}

// ---- The model document: facts that trace-stream needs beyond the compiled model ----

[[nodiscard]] inline const Value &at(const Value &object, std::string_view key) {
  const Value *found = object.find(key);
  if (found == nullptr)
    throw ModelError("optical model: missing '" + std::string(key) + "'");
  return *found;
}

[[nodiscard]] inline double number(const Value &value) {
  if (value.kind() == Value::Kind::boolean)
    return std::get<bool>(value.data) ? 1.0 : 0.0;
  if (value.kind() != Value::Kind::number)
    throw ModelError("optical model: expected a number");
  return std::get<double>(value.data);
}

[[nodiscard]] inline const Value::Array &items(const Value &value) {
  if (value.kind() != Value::Kind::array)
    throw ModelError("optical model: expected an array");
  return std::get<Value::Array>(value.data);
}

[[nodiscard]] inline Vec vec3(const Value &value) {
  const auto &list = items(value);
  if (list.size() != 3)
    throw ModelError("optical model: expected three components");
  return {number(list[0]), number(list[1]), number(list[2])};
}

[[nodiscard]] inline std::string text(const Value &value) {
  if (value.kind() != Value::Kind::string)
    throw ModelError("optical model: expected a string");
  return std::get<std::string>(value.data);
}

[[nodiscard]] inline bool truthy(const Value *value) {
  if (value == nullptr)
    return false;
  switch (value->kind()) {
  case Value::Kind::null:
    return false;
  case Value::Kind::boolean:
    return std::get<bool>(value->data);
  case Value::Kind::number:
    return std::get<double>(value->data) != 0.0;
  case Value::Kind::string:
    return !std::get<std::string>(value->data).empty();
  case Value::Kind::array:
    return !std::get<Value::Array>(value->data).empty();
  case Value::Kind::object:
    return !std::get<Value::Object>(value->data).names.empty();
  }
  return false;
}

[[nodiscard]] inline bool has_leg(const Value &obscurer, std::string_view leg) {
  const Value *legs = obscurer.find("legs");
  if (legs == nullptr)
    return false;
  for (const auto &item : items(*legs))
    if (item.kind() == Value::Kind::string && std::get<std::string>(item.data) == leg)
      return true;
  return false;
}

// Radius [m] of the circle that encloses a plate aperture, as in the native tracer.
[[nodiscard]] inline double circumradius(const std::string &aperture, double diameter_m) {
  const double half_width = 0.5 * diameter_m;
  if (aperture == "square")
    return half_width * std::sqrt(2.0);
  if (aperture.rfind("hexagon", 0) == 0)
    return half_width * 2.0 / std::sqrt(3.0);
  return half_width;
}

[[nodiscard]] inline std::array<double, 2> end_radii(const Value &obscurer) {
  if (obscurer.find("start_radius_m") != nullptr) {
    const double thickness = number(at(obscurer, "thickness_m"));
    return {number(at(obscurer, "start_radius_m")) + thickness,
            number(at(obscurer, "end_radius_m")) + thickness};
  }
  const double radius = 0.5 * number(at(obscurer, "diameter_m"));
  return {radius, radius};
}

[[nodiscard]] inline Vec unit(const Vec &v) {
  const double norm = std::sqrt(dot(v, v));
  return {v[0] / norm, v[1] / norm, v[2] / norm};
}

[[nodiscard]] inline Vec minus(const Vec &a, const Vec &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

// Highest z [m] of one obscurer record.
[[nodiscard]] inline double obscurer_top_z(const Value &obscurer) {
  if (const Value *vertices = obscurer.find("vertices_m")) {
    double top = -std::numeric_limits<double>::infinity();
    for (const auto &vertex : items(*vertices))
      top = std::max(top, vec3(vertex)[2]);
    return top;
  }
  if (obscurer.find("centre_m") != nullptr) {
    const Vec normal = unit(vec3(at(obscurer, "normal")));
    const double radius =
        circumradius(text(at(obscurer, "aperture")), number(at(obscurer, "diameter_m")));
    return vec3(at(obscurer, "centre_m"))[2] +
           radius * std::sqrt(std::max(0.0, 1.0 - std::pow(normal[2], 2.0)));
  }
  const Vec start = vec3(at(obscurer, "start_m")), end = vec3(at(obscurer, "end_m"));
  const Vec axis = unit(minus(end, start));
  const double sine = std::sqrt(std::max(0.0, 1.0 - std::pow(axis[2], 2.0)));
  const auto radii = end_radii(obscurer);
  return std::max(start[2] + radii[0] * sine, end[2] + radii[1] * sine);
}

// Maximum dot product with a direction over an obscurer's geometry.
[[nodiscard]] inline double obscurer_support(const Value &obscurer, const Vec &direction) {
  if (const Value *vertices = obscurer.find("vertices_m")) {
    double best = -std::numeric_limits<double>::infinity();
    for (const auto &vertex : items(*vertices))
      best = std::max(best, dot(vec3(vertex), direction));
    return best;
  }
  const auto radial_of = [&](const Vec &axis) {
    return std::sqrt(
        std::max(0.0, dot(direction, direction) - std::pow(dot(direction, axis), 2.0)));
  };
  if (obscurer.find("centre_m") != nullptr) {
    const Vec normal = unit(vec3(at(obscurer, "normal")));
    const double radius =
        circumradius(text(at(obscurer, "aperture")), number(at(obscurer, "diameter_m")));
    return dot(vec3(at(obscurer, "centre_m")), direction) + radius * radial_of(normal);
  }
  const Vec start = vec3(at(obscurer, "start_m")), end = vec3(at(obscurer, "end_m"));
  const double radial = radial_of(unit(minus(end, start)));
  const auto radii = end_radii(obscurer);
  return std::max(dot(start, direction) + radii[0] * radial,
                  dot(end, direction) + radii[1] * radial);
}

// Conservative bound on the dot product over an asphere's aperture (defined below).
[[nodiscard]] inline double asphere_support(const Value &asphere, const Vec &direction);

// Conservative bound on the dot product with a direction over a mirror's surface.
[[nodiscard]] inline double mirror_support(const Value &mirror, const Vec &direction) {
  if (const Value *facets = mirror.find("facets")) {
    double best = -std::numeric_limits<double>::infinity();
    for (const auto &facet : items(*facets)) {
      const Value *curvature = facet.find("radius_of_curvature_m");
      const double radius = curvature == nullptr ? 0.0 : number(*curvature);
      const double reach = circumradius(text(at(facet, "shape")), number(at(facet, "diameter_m")));
      const double sag =
          radius != 0.0 ? radius - std::sqrt(std::max(0.0, radius * radius - reach * reach)) : 0.0;
      best = std::max(best, dot(vec3(at(facet, "centre_m")), direction) + reach + sag);
    }
    return best;
  }
  return asphere_support(at(mirror, "asphere"), direction);
}

// Conservative bound on the dot product over an asphere's aperture (a mirror's or the focal
// surface's `asphere` record).
[[nodiscard]] inline double asphere_support(const Value &asphere, const Vec &direction) {
  const double inner = number(at(asphere, "inner_radius_m"));
  const double outer = number(at(asphere, "outer_radius_m"));
  const double scale = number(at(asphere, "radial_scale_m"));
  const double vertex = number(at(asphere, "vertex_z_m"));
  const auto &coefficients = items(at(asphere, "coefficients_m"));
  // numpy.linspace(inner, outer, kAsphereSamples).
  const double step = (outer - inner) / static_cast<double>(kAsphereSamples - 1);
  std::vector<double> radii(kAsphereSamples), z(kAsphereSamples, vertex);
  for (std::size_t k = 0; k < kAsphereSamples; ++k)
    radii[k] = static_cast<double>(k) * step + inner;
  radii.back() = outer;
  double curvature_bound = 0.0;
  double magnitude = std::abs(vertex);
  const double outer_scaled = outer / scale;
  for (std::size_t i = 0; i < coefficients.size(); ++i) {
    const double coefficient = number(coefficients[i]);
    const double power = static_cast<double>(2 * i);
    for (std::size_t k = 0; k < kAsphereSamples; ++k)
      z[k] += coefficient * std::pow(radii[k] / scale, power);
    magnitude += std::abs(coefficient) * std::pow(outer_scaled, power);
    if (i != 0)
      curvature_bound += std::abs(coefficient) * power * (power - 1.0) *
                         std::pow(outer_scaled, power - 2.0) / std::pow(scale, 2.0);
  }
  const double radial = std::hypot(direction[0], direction[1]);
  const double dz = std::abs(direction[2]);
  const double rounding = std::numeric_limits<double>::epsilon() *
                          static_cast<double>(16 + 4 * coefficients.size()) *
                          (radial * outer + dz * magnitude);
  double best = -std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < kAsphereSamples; ++k)
    best = std::max(best, radial * radii[k] + direction[2] * z[k]);
  return best + (dz * curvature_bound * std::pow(step, 2.0) / 8.0 + rounding);
}

// An optical model with the document facts that trace-stream reads (OpticalModel in Python).
struct StreamModel {
  OpticalModel model;
  Value document;
  std::string sha256;
  bool has_structure{};
  bool has_camera_window{};                   // optics.camera_window: the window is traced
  std::optional<double> transmission_on_axis; // optics.telescope_transmission.on_axis

  [[nodiscard]] const Value &optics() const { return at(document, "optics"); }

  [[nodiscard]] Vec axes_offset() const {
    const Value *offset = optics().find("axes_offset_m");
    return offset == nullptr ? Vec{0.0, 0.0, 0.0} : vec3(*offset);
  }

  [[nodiscard]] std::vector<const Value *> incoming_obscurers() const {
    std::vector<const Value *> out;
    for (const auto &obscurer : items(at(optics(), "obscurers")))
      if (has_leg(obscurer, "incoming"))
        out.push_back(&obscurer);
    return out;
  }

  // The primary, and the secondary if it blocks incoming light: (role, support).
  [[nodiscard]] std::vector<std::pair<std::string, double>>
  mirrors_support(const Vec &direction) const {
    std::vector<std::string> roles = {"primary"};
    if (truthy(optics().find("secondary_blocks_incoming")))
      roles.emplace_back("secondary");
    std::vector<std::pair<std::string, double>> out;
    for (const auto &role : roles) {
      const Value *mirror = optics().find(role);
      if (mirror != nullptr && mirror->kind() != Value::Kind::null)
        out.emplace_back(role, mirror_support(*mirror, direction));
    }
    return out;
  }

  // Highest incoming obscurer top [m] and its name.
  [[nodiscard]] std::optional<std::pair<double, std::string>> incoming_top() const {
    std::optional<std::pair<double, std::string>> top;
    for (const Value *obscurer : incoming_obscurers()) {
      std::pair<double, std::string> item{obscurer_top_z(*obscurer), text(at(*obscurer, "name"))};
      if (!top || item > *top)
        top = std::move(item);
    }
    return top;
  }

  [[nodiscard]] std::optional<double> max_mirror_z() const {
    std::optional<double> top;
    for (const auto &[role, support] : mirrors_support({0.0, 0.0, 1.0}))
      top = top ? std::max(*top, support) : support;
    return top;
  }
};

// Largest distance [m] from `centre` to a body, bounded through its support function (an upper
// bound of the dot product with any direction): the body lies in the box spanned by its supports
// along the three axes, so its distance from the centre is at most the box's diagonal reach.
template <class Support>
[[nodiscard]] inline double support_box_distance(const Support &support, const Vec &centre) {
  double sum = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    Vec axis{0.0, 0.0, 0.0};
    axis[i] = 1.0;
    const Vec flipped{-axis[0], -axis[1], -axis[2]};
    const double reach = std::max(support(axis) - centre[i], support(flipped) + centre[i]);
    sum += reach * reach;
  }
  return std::sqrt(sum);
}

// Name and distance [m] from the sphere centre (in the traced frame, -axes_offset) of every
// mirror that exists, every obscurer (any leg) and the focal surface. The distance of a mirror,
// an obscurer or an asphere focal surface is an upper bound (support box). The focal plane has no
// lateral extent in the model: its entry is the distance of its height from the centre, so the
// check requires the plane to reach the sphere.
[[nodiscard]] inline std::vector<std::pair<std::string, double>>
enclosure_distances(const StreamModel &model) {
  const Vec centre = minus({0.0, 0.0, 0.0}, model.axes_offset());
  std::vector<std::pair<std::string, double>> out;
  for (const std::string role : {"primary", "secondary"}) {
    const Value *mirror = model.optics().find(role);
    if (mirror != nullptr && mirror->kind() != Value::Kind::null)
      out.emplace_back(role, support_box_distance(
                                 [&](const Vec &d) { return mirror_support(*mirror, d); }, centre));
  }
  for (const auto &obscurer : items(at(model.optics(), "obscurers")))
    out.emplace_back(
        text(at(obscurer, "name")),
        support_box_distance([&](const Vec &d) { return obscurer_support(obscurer, d); }, centre));
  const Value &focal = at(model.optics(), "focal_surface");
  const Value *asphere = focal.find("asphere");
  if (asphere != nullptr && asphere->kind() != Value::Kind::null)
    out.emplace_back(
        "focal surface",
        support_box_distance([&](const Vec &d) { return asphere_support(*asphere, d); }, centre));
  else
    out.emplace_back("focal surface (plane)", std::abs(number(at(focal, "plane_z_m")) - centre[2]));
  return out;
}

// Loads an optical-model file's text; the sha256 entry must match the content.
[[nodiscard]] inline std::shared_ptr<const StreamModel> load_stream_model(std::string_view text) {
  auto out = std::make_shared<StreamModel>();
  try {
    out->document = json::Parser::parse(text);
  } catch (const std::runtime_error &error) {
    throw ModelError(error.what());
  }
  if (out->document.kind() != Value::Kind::object)
    throw ModelError("optical model: expected a JSON object");
  const Value *sha = out->document.find("sha256");
  if (sha == nullptr || sha->kind() != Value::Kind::string ||
      std::get<std::string>(sha->data) != model_hash(text))
    throw ModelError("optical model content does not match its sha256");
  out->sha256 = std::get<std::string>(sha->data);
  out->model = load_model(text);
  out->has_structure = truthy(&at(at(out->document, "provenance"), "structure"));
  out->has_camera_window = out->model.camera_window.has_value();
  if (const Value *transmission = out->optics().find("telescope_transmission");
      transmission != nullptr && transmission->kind() != Value::Kind::null)
    out->transmission_on_axis = number(at(*transmission, "on_axis"));
  return out;
}

// Rotation R with p_telescope = R p_ground for a telescope at zenith theta, azimuth phi
// (row-major; see ground_to_telescope_rotation in Python).
[[nodiscard]] inline std::array<double, 9> ground_to_telescope_rotation(double theta_deg,
                                                                        double phi_deg) {
  const double theta = theta_deg * (kPi / 180.0), phi = phi_deg * (kPi / 180.0);
  // M = step(alpha = phi; 0, 1) @ step(beta = -theta; 2, 0) for the row-vector rule; R = M^T.
  const auto step = [](double angle, int i, int j) {
    std::array<double, 9> m = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    m[3 * i + i] = m[3 * j + j] = std::cos(angle);
    m[3 * i + j] = std::sin(angle);
    m[3 * j + i] = -std::sin(angle);
    return m;
  };
  const auto a = step(phi, 0, 1), b = step(-theta, 2, 0);
  std::array<double, 9> m{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      m[3 * i + j] = a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j];
  std::array<double, 9> r{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r[3 * i + j] = m[3 * j + i];
  return r;
}

[[nodiscard]] inline Vec column(const std::array<double, 9> &r, int j) {
  return {r[j], r[3 + j], r[6 + j]};
}

// Height [m] of the ground-frame input plane above the telescope sphere centre.
[[nodiscard]] inline double minimum_input_plane_height_m(const StreamModel &model, double theta_deg,
                                                         double phi_deg) {
  if (!(theta_deg >= 0.0 && theta_deg < 90.0))
    throw ModelError("telescope zenith angle " + std::string(format_g(theta_deg)) +
                     " deg must lie in [0, 90)");
  const Vec normal = column(ground_to_telescope_rotation(theta_deg, phi_deg), 2);
  std::optional<double> support;
  for (const Value *obscurer : model.incoming_obscurers()) {
    const double s = obscurer_support(*obscurer, normal);
    support = support ? std::max(*support, s) : s;
  }
  for (const auto &[role, s] : model.mirrors_support(normal))
    support = support ? std::max(*support, s) : s;
  return std::max(0.0, support.value_or(0.0) + dot(normal, model.axes_offset())) +
         kInputPlaneMarginM;
}

// Every obscurer (any leg) and the primary and secondary mirrors, with the lowest point of each
// along `normal`: min n.p over the optic, from the support in direction -normal. Photons from
// below meet all of them.
[[nodiscard]] inline std::vector<std::pair<std::string, double>>
optics_bottoms(const StreamModel &model, const Vec &normal) {
  const Vec flipped = {-normal[0], -normal[1], -normal[2]};
  std::vector<std::pair<std::string, double>> out;
  for (const auto &obscurer : items(at(model.optics(), "obscurers")))
    out.emplace_back(text(at(obscurer, "name")), -obscurer_support(obscurer, flipped));
  for (const std::string role : {"primary", "secondary"}) {
    const Value *mirror = model.optics().find(role);
    if (mirror != nullptr && mirror->kind() != Value::Kind::null)
      out.emplace_back(role, -mirror_support(*mirror, flipped));
  }
  return out;
}

// Input plane of photons from below, in the input frame (the same coordinate as the upper
// plane): the plane lies a margin below the lowest optic along `normal`, with the traced origin
// offset `traced_offset` (offset - axes_offset) as for the upper plane.
[[nodiscard]] inline double lower_input_plane_m(const StreamModel &model, const Vec &normal,
                                                const Vec &traced_offset) {
  std::optional<double> bottom;
  for (const auto &[name, lowest] : optics_bottoms(model, normal))
    bottom = bottom ? std::min(*bottom, lowest) : lowest;
  return bottom.value_or(0.0) - kInputPlaneMarginM - dot(normal, traced_offset);
}

// The input-plane checks of OpticalModel.trace_bunches (Python): the plane must lie above the
// incoming obscurers and above the mirrors that photons meet first.
inline void check_input_plane(const StreamModel &model, const std::array<double, 9> &rotation,
                              const Vec &offset, double plane_z) {
  const Vec axes = model.axes_offset();
  const Vec traced_offset = minus(offset, axes);
  const Vec normal = column(rotation, 2);
  const Vec point = {traced_offset[0] + normal[0] * plane_z * 1.0,
                     traced_offset[1] + normal[1] * plane_z * 1.0,
                     traced_offset[2] + normal[2] * plane_z * 1.0};
  const double value = dot(normal, point);
  if (const auto top = model.incoming_top()) {
    std::vector<std::tuple<double, double, std::string>> checks;
    if (normal[2] > 0.0) {
      for (const Value *obscurer : model.incoming_obscurers())
        checks.emplace_back(value / normal[2], obscurer_support(*obscurer, normal) / normal[2],
                            text(at(*obscurer, "name")));
    } else {
      checks.emplace_back(-std::numeric_limits<double>::infinity(), top->first, top->second);
    }
    for (const auto &[height, obscurer_top, name] : checks)
      if (height < obscurer_top)
        throw ModelError("input plane below incoming obscurer '" + name + "': the plane is at " +
                         format_3f(height) + " m along the telescope axis, but the obscurer " +
                         "reaches " + format_3f(obscurer_top) +
                         " m; raise the input plane above it");
  }
  for (const auto &[role, mirror_top] : model.mirrors_support(normal))
    if (value < mirror_top)
      throw ModelError("input plane below the " + role + " mirror: the plane is at n.p = " +
                       format_3f(value) + " m but the mirror reaches " + format_3f(mirror_top) +
                       " m along the input normal; raise the input plane above it");
}

// The check of the plane for photons from below: it must lie below every optic (any leg).
inline void check_lower_input_plane(const StreamModel &model, const std::array<double, 9> &rotation,
                                    const Vec &offset, double plane_z) {
  const Vec traced_offset = minus(offset, model.axes_offset());
  const Vec normal = column(rotation, 2);
  const Vec point = {traced_offset[0] + normal[0] * plane_z, traced_offset[1] + normal[1] * plane_z,
                     traced_offset[2] + normal[2] * plane_z};
  const double value = dot(normal, point);
  for (const auto &[name, bottom] : optics_bottoms(model, normal))
    if (value > bottom)
      throw ModelError("input plane above the " + name +
                       " for photons from below: the plane is "
                       "at n.p = " +
                       format_3f(value) + " m but the optic reaches down to " + format_3f(bottom) +
                       " m along the input normal; lower the input plane "
                       "below it");
}

// ---- Efficiency status ----

[[nodiscard]] inline std::vector<std::string> names_of(const st::Metadata &metadata) {
  std::vector<std::string> names;
  const auto found = metadata.find("factors");
  if (found == metadata.end())
    return names;
  std::string text;
  if (const auto *s = std::get_if<std::string>(&found->second))
    text = *s;
  else if (const auto *i = std::get_if<std::int64_t>(&found->second))
    text = std::to_string(*i);
  std::size_t start = 0;
  while (start <= text.size()) {
    const std::size_t comma = std::min(text.find(',', start), text.size());
    std::string item = text.substr(start, comma - start);
    const auto first = item.find_first_not_of(" \t\n\r\f\v");
    if (first != std::string::npos) {
      const auto last = item.find_last_not_of(" \t\n\r\f\v");
      names.push_back(item.substr(first, last - first + 1));
    }
    start = comma + 1;
  }
  return names;
}

[[nodiscard]] inline std::string join(const std::vector<std::string> &names) {
  std::string out;
  for (std::size_t i = 0; i < names.size(); ++i)
    out += (i ? "," : "") + names[i];
  return out;
}

[[nodiscard]] inline bool contains(const std::vector<std::string> &names, std::string_view name) {
  return std::find(names.begin(), names.end(), name) != names.end();
}

// One efficiency status entry; for `envelope`, its efficiency bound table.
struct Factor {
  std::string state;
  std::optional<st::Value> wavelength_nm;
  std::optional<st::Value> value;
};

[[nodiscard]] inline std::string value_text(const st::Value &value) {
  if (const auto *s = std::get_if<std::string>(&value))
    return *s;
  if (const auto *i = std::get_if<std::int64_t>(&value))
    return std::to_string(*i);
  if (const auto *d = std::get_if<double>(&value))
    return hash_detail::python_float_repr(*d);
  return "[...]";
}

[[nodiscard]] inline std::map<std::string, Factor> factors_of(const st::Metadata &metadata) {
  std::map<std::string, Factor> factors;
  for (const auto &name : names_of(metadata)) {
    Factor factor;
    const auto state = metadata.find(name + ".state");
    factor.state = state == metadata.end() ? "pending" : value_text(state->second);
    if (factor.state != "applied" && factor.state != "envelope" && factor.state != "pending")
      throw ModelError("efficiency_status factor '" + name + "': unknown state '" + factor.state +
                       "'");
    if (const auto w = metadata.find(name + ".wavelength_nm"); w != metadata.end())
      factor.wavelength_nm = w->second;
    if (const auto v = metadata.find(name + ".value"); v != metadata.end())
      factor.value = v->second;
    factors.insert_or_assign(name, std::move(factor));
  }
  return factors;
}

[[nodiscard]] inline std::vector<double> as_array(const st::Value &value) {
  if (const auto *v = std::get_if<std::vector<double>>(&value))
    return *v;
  if (const auto *d = std::get_if<double>(&value))
    return {*d};
  if (const auto *i = std::get_if<std::int64_t>(&value))
    return {static_cast<double>(*i)};
  double parsed = 0.0;
  if (obdeect::stream::parse::finite_number(std::get<std::string>(value), parsed))
    return {parsed};
  throw ModelError("could not convert string to float: '" + std::get<std::string>(value) + "'");
}

struct Envelope {
  std::vector<double> nodes, values;
};

[[nodiscard]] inline Envelope envelope_table(const std::string &name, const Factor &factor) {
  if (!factor.wavelength_nm || !factor.value)
    throw ModelError("efficiency status factor '" + name + "' is 'envelope' but has no " + name +
                     ".wavelength_nm and " + name + ".value");
  Envelope table{as_array(*factor.wavelength_nm), as_array(*factor.value)};
  bool valid = !table.nodes.empty() && table.nodes.size() == table.values.size();
  for (std::size_t i = 1; valid && i < table.nodes.size(); ++i)
    valid = table.nodes[i] - table.nodes[i - 1] > 0.0;
  for (std::size_t i = 0; valid && i < table.values.size(); ++i)
    valid = table.values[i] > 0.0;
  if (!valid)
    throw ModelError("efficiency_status factor '" + name +
                     "': the envelope needs increasing wavelengths and positive values of the "
                     "same length");
  return table;
}

// ---- The tracer ----

enum class TransmissionPart { full, angular };

// Optical model of a telescope, with the part of its transmission that the trace applies.
struct Optics {
  std::shared_ptr<const StreamModel> model;
  TransmissionPart part{TransmissionPart::full};
  double on_axis{};
  double telescope_plane_m{};       // Input plane of telescope-frame photons.
  double telescope_lower_plane_m{}; // Input plane of telescope-frame photons from below.
  // Sphere runs: name and distance [m] from the sphere centre of each optical element (see
  // enclosure_distances); the sphere of a telescope must reach every one.
  std::vector<std::pair<std::string, double>> enclosure;
  std::string provenance;

  [[nodiscard]] bool normalise() const { return part == TransmissionPart::angular; }

  [[nodiscard]] std::vector<std::string> own_factors() const {
    std::vector<std::string> factors{"mirror_reflectivity"};
    if (model->has_structure || model->transmission_on_axis)
      factors.emplace_back("telescope_transmission");
    // A traced camera window (Fresnel and bulk loss); sim_telarray never applies it.
    if (model->has_camera_window)
      factors.emplace_back("camera_window");
    return factors;
  }
};

[[nodiscard]] inline std::shared_ptr<const Optics>
make_optics(std::shared_ptr<const StreamModel> model, TransmissionPart part,
            const std::string &version) {
  auto optics = std::make_shared<Optics>();
  optics->part = part;
  if (part == TransmissionPart::angular) {
    if (!model->transmission_on_axis)
      throw ModelError("the angular telescope transmission needs a model with "
                       "telescope_transmission; this model has none");
    optics->on_axis = *model->transmission_on_axis;
    if (!(optics->on_axis > 0.0))
      throw ModelError("the angular telescope transmission needs a positive on-axis value");
  }
  std::optional<double> top;
  if (const auto obscurer = model->incoming_top())
    top = obscurer->first;
  if (const auto mirror = model->max_mirror_z())
    top = top ? std::max(*top, *mirror) : *mirror;
  optics->telescope_plane_m = top.value_or(0.0) + kInputPlaneMarginM;
  optics->telescope_lower_plane_m = lower_input_plane_m(*model, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
  optics->enclosure = enclosure_distances(*model);
  optics->provenance = "obdeect " + version + " model " + model->sha256;
  optics->model = std::move(model);
  return optics;
}

// The efficiency status record with the factors obdeect applies marked applied.
[[nodiscard]] inline st::Metadata with_applied(const st::Metadata &metadata, const Optics &optics) {
  const auto all_names = names_of(metadata);
  if (optics.normalise() && contains(all_names, kConstantFactor))
    throw ModelError(std::string("the efficiency status already names ") + kConstantFactor);
  std::vector<std::string> names;
  const auto own = optics.own_factors();
  for (const auto &name : all_names)
    if (contains(own, name))
      names.push_back(name);
  if (names.empty())
    return metadata;
  st::Metadata updated = metadata;
  for (const auto &name : names) {
    updated[name + ".state"] = std::string("applied");
    updated[name + ".source"] = optics.provenance;
  }
  if (contains(names, "telescope_transmission")) {
    updated["telescope_transmission.part"] = std::string(optics.normalise() ? "angular" : "full");
    if (optics.normalise()) {
      updated[std::string(kConstantFactor) + ".state"] = std::string("pending");
      updated[std::string(kConstantFactor) + ".value"] = optics.on_axis;
      auto factors = names_of(updated);
      factors.emplace_back(kConstantFactor);
      updated["factors"] = join(factors);
    }
  }
  return updated;
}

// The efficiency status that marks the factors obdeect applies as applied, merged into an
// incoming one (each factor named once).
[[nodiscard]] inline st::Metadata own_status(const Optics &optics,
                                             const st::Metadata &metadata = {}) {
  st::Metadata merged = metadata;
  auto names = names_of(metadata);
  for (const auto &name : optics.own_factors())
    if (!contains(names, name))
      names.push_back(name);
  merged["factors"] = join(names);
  return with_applied(merged, optics);
}

[[nodiscard]] inline std::int64_t int_key(const st::Metadata &metadata, const std::string &key,
                                          const std::string &record) {
  const auto found = metadata.find(key);
  if (found == metadata.end())
    throw ModelError(record + " record without " + key);
  if (const auto *i = std::get_if<std::int64_t>(&found->second))
    return *i;
  if (const auto *d = std::get_if<double>(&found->second))
    return static_cast<std::int64_t>(*d);
  if (const auto *s = std::get_if<std::string>(&found->second)) {
    std::int64_t value = 0;
    if (obdeect::stream::parse::integer(*s, value))
      return value;
  }
  throw ModelError(record + " record: " + key + " is not an integer");
}

[[nodiscard]] inline double float_key(const st::Value &value, const std::string &key) {
  if (const auto *d = std::get_if<double>(&value))
    return *d;
  if (const auto *i = std::get_if<std::int64_t>(&value))
    return static_cast<double>(*i);
  if (const auto *s = std::get_if<std::string>(&value)) {
    double parsed = 0.0;
    if (obdeect::stream::parse::finite_number(*s, parsed))
      return parsed;
  }
  throw ModelError("run record: " + key + " is not a number");
}

struct Config {
  bool roulette{true};
  std::vector<Envelope> envelopes;
};

// The rotation to the telescope frame, the offset of the traced frame and the input planes of
// the photons of one telescope pointing (or of the telescope frame): upper for photons from
// above, lower for photons from below. The plane checks depend on these values only, so they run
// once (at the first photons that need them).
struct InputPlanes {
  std::array<double, 9> rotation{1, 0, 0, 0, 1, 0, 0, 0, 1};
  Vec offset{0.0, 0.0, 0.0};
  double upper{};
  double lower{};
  bool upper_checked{};
  bool lower_checked{};
};

// Ground-frame planes of a telescope pointing (zenith, azimuth) [deg].
[[nodiscard]] inline InputPlanes ground_input_planes(const StreamModel &model, double zenith_deg,
                                                     double azimuth_deg) {
  InputPlanes planes;
  planes.rotation = ground_to_telescope_rotation(zenith_deg, azimuth_deg);
  // The traced origin is -axes_offset (the ground origin has no offset).
  planes.upper = minimum_input_plane_height_m(model, zenith_deg, azimuth_deg);
  planes.lower = lower_input_plane_m(model, column(planes.rotation, 2),
                                     minus({0.0, 0.0, 0.0}, model.axes_offset()));
  return planes;
}

// Telescope-frame planes: no tilt. The traced frame is that of trace(), so the axes offset that
// it subtracts is cancelled, and a photon at p is traced at p.
[[nodiscard]] inline InputPlanes telescope_input_planes(const Optics &optics) {
  InputPlanes planes;
  planes.offset = optics.model->axes_offset();
  planes.upper = optics.telescope_plane_m;
  planes.lower = optics.telescope_lower_plane_m;
  return planes;
}

// Per-telescope state of the block that is being read.
struct Telescope {
  std::uint32_t id{};
  bool kept{};
  std::map<std::string, Factor> factors;
  std::optional<Config> config;
  bool photons_seen{};
  bool status_written{};
  std::uint64_t out_index{};
  double weight{};
};

} // namespace stream_trace_detail

// Settings of one trace-stream run (see `obdeect-ray-tracing-trace-stream`).
struct StreamTraceSettings {
  std::uint64_t seed{};
  unsigned threads{1};
  std::optional<std::set<std::uint32_t>> telescopes; // None: all (or all with a model)
  std::optional<std::pair<double, double>> pointing; // (zenith, azimuth) [deg]
  stream_trace_detail::TransmissionPart part{stream_trace_detail::TransmissionPart::full};
  std::string version; // In the producer and provenance strings.
};

struct StreamTraceCounts {
  std::uint64_t events{}, photons_in{}, photons_out{};
};

// One optical model of a run, as the text of its JSON file: the model of every telescope (no
// telescope id), or the model of one telescope id.
struct ModelSource {
  std::optional<std::uint32_t> telescope;
  std::string text;
};

// A telescope id given by the user: sim_telarray numbering from 1; the random key holds 16 bits.
[[nodiscard]] inline std::uint32_t checked_telescope_id(std::int64_t value) {
  if (value < 1 || value > std::numeric_limits<std::uint16_t>::max())
    throw ModelError("telescope id " + std::to_string(value) + " outside 1.." +
                     std::to_string(std::numeric_limits<std::uint16_t>::max()));
  return static_cast<std::uint32_t>(value);
}

// The optics of a run and the telescopes it traces, from its models: one model for every
// telescope, or the model of each telescope id. A model text is loaded once.
struct RunOptics {
  std::shared_ptr<const stream_trace_detail::Optics> shared;
  std::map<std::uint32_t, std::shared_ptr<const stream_trace_detail::Optics>> per_telescope;
  std::optional<std::set<std::uint32_t>> telescopes; // None: every telescope
};

[[nodiscard]] inline RunOptics make_run_optics(const std::vector<ModelSource> &models,
                                               std::optional<std::set<std::uint32_t>> telescopes,
                                               stream_trace_detail::TransmissionPart part,
                                               const std::string &version) {
  RunOptics run;
  if (models.empty())
    throw ModelError("give at least one optical model");
  if (models.size() == 1 && !models.front().telescope) {
    run.shared = stream_trace_detail::make_optics(
        stream_trace_detail::load_stream_model(models.front().text), part, version);
    run.telescopes = std::move(telescopes);
    return run;
  }
  std::map<std::string, std::shared_ptr<const stream_trace_detail::Optics>> loaded;
  for (const auto &model : models) {
    if (!model.telescope)
      throw ModelError("give one optical model for every telescope, or one for each telescope id");
    if (run.per_telescope.count(*model.telescope))
      throw ModelError("a telescope id is given in more than one optical model");
    auto &optics = loaded[model.text];
    if (!optics)
      optics = stream_trace_detail::make_optics(stream_trace_detail::load_stream_model(model.text),
                                                part, version);
    run.per_telescope[*model.telescope] = optics;
  }
  std::set<std::uint32_t> selected;
  if (telescopes) {
    selected = *telescopes;
  } else {
    for (const auto &[id, optics] : run.per_telescope)
      selected.insert(id);
  }
  std::string missing;
  for (const auto id : selected)
    if (!run.per_telescope.count(id))
      missing += (missing.empty() ? "" : ", ") + std::to_string(id);
  if (!missing.empty())
    throw ModelError("telescopes without an optical model: [" + missing + "]");
  for (auto it = run.per_telescope.begin(); it != run.per_telescope.end();)
    it = selected.count(it->first) ? std::next(it) : run.per_telescope.erase(it);
  run.telescopes = std::move(selected);
  return run;
}

// True if two paths name the same file (through links): an input must not be overwritten.
[[nodiscard]] inline bool same_file(const std::string &input, const std::string &output) {
  std::error_code error;
  const bool same = std::filesystem::equivalent(input, output, error);
  return !error && same;
}

// Reads the records of one stream and writes the traced records.
class StreamTracer {
public:
  using Optics = stream_trace_detail::Optics;

  // `shared`: one model for every telescope; `per_telescope`: the model of each telescope id.
  StreamTracer(std::shared_ptr<const Optics> shared,
               std::map<std::uint32_t, std::shared_ptr<const Optics>> per_telescope,
               obdeect::stream::Writer &writer, const StreamTraceSettings &settings)
      : shared_(std::move(shared)), per_telescope_(std::move(per_telescope)), writer_(writer),
        settings_(settings), pool_(settings.threads) {}

  void handle(const obdeect::stream::Record &record) {
    using obdeect::stream::RecordType;
    const RecordType kind = record.type;
    if (run_status_pending_ && (kind == RecordType::event || kind == RecordType::run_end)) {
      writer_.write_metadata(RecordType::efficiency_status,
                             stream_trace_detail::own_status(*shared_));
      run_status_pending_ = false;
    }
    switch (kind) {
    case RecordType::run:
      run(record.metadata);
      break;
    case RecordType::efficiency_status:
      efficiency_status(record.metadata);
      break;
    case RecordType::event:
      event(record.metadata);
      break;
    case RecordType::telescope:
      telescope(record.metadata);
      break;
    case RecordType::photons:
      photons(record.chunk);
      break;
    case RecordType::telescope_end:
      telescope_end(record.metadata);
      break;
    case RecordType::blob:
      writer_.write_blob(record.blob_eventio_type, record.blob_eventio_id, record.blob);
      break;
    case RecordType::event_end:
    case RecordType::run_end:
      writer_.write_metadata(kind, record.metadata);
      break;
    case RecordType::stream:
    case RecordType::stream_end:
      break;
    default:
      throw ModelError(std::string("unexpected photon-stream record '") +
                       obdeect::stream::record_name(kind) + "'");
    }
  }

  [[nodiscard]] const StreamTraceCounts &counts() const { return counts_; }

private:
  using Metadata = obdeect::stream::Metadata;
  using Telescope = stream_trace_detail::Telescope;
  using InputPlanes = stream_trace_detail::InputPlanes;
  // Input planes by optical model and pointing (none: the telescope frame).
  using GeometryKey = std::tuple<const Optics *, std::optional<std::pair<double, double>>>;

  [[nodiscard]] bool per_telescope() const { return !shared_; }

  [[nodiscard]] const Optics &optics(std::uint32_t telescope_id) const {
    if (shared_)
      return *shared_;
    return *per_telescope_.at(telescope_id);
  }

  void run(const Metadata &metadata) {
    using namespace stream_trace_detail;
    run_factors_.clear();
    law_.reset();
    if (metadata.count("wavelength_min_nm") && metadata.count("wavelength_max_nm")) {
      const double lo = float_key(metadata.at("wavelength_min_nm"), "wavelength_min_nm");
      const double hi = float_key(metadata.at("wavelength_max_nm"), "wavelength_max_nm");
      // One wavelength: the 1/lambda^2 law over [lo, lo + 1) truncates every draw to floor(lo).
      law_ = InverseSquareLaw{lo, hi == lo ? lo + 1.0 : hi};
    }
    const bool zenith = metadata.count("telescope_pointing_zenith_rad") != 0;
    const bool azimuth = metadata.count("telescope_pointing_azimuth_rad") != 0;
    run_pointing_.reset();
    if (zenith != azimuth)
      throw ModelError("run record: give both telescope_pointing_zenith_rad and "
                       "telescope_pointing_azimuth_rad or neither");
    if (zenith) {
      const double radians_to_degrees = 180.0 / kPi;
      run_pointing_ = std::pair{
          float_key(metadata.at("telescope_pointing_zenith_rad"), "telescope_pointing_zenith_rad") *
              radians_to_degrees,
          float_key(metadata.at("telescope_pointing_azimuth_rad"),
                    "telescope_pointing_azimuth_rad") *
              radians_to_degrees};
    }
    // Ground representation (format 1.2): sphere runs have G rows at their own 3D points.
    try {
      sphere_ = st::ground_representation(metadata) == st::GroundRepresentation::sphere;
    } catch (const std::runtime_error &error) {
      throw ModelError(error.what());
    }
    radii_.clear();
    if (sphere_) {
      const auto found = metadata.find("telescope_r_m");
      if (found != metadata.end()) {
        const auto *radii = std::get_if<std::vector<double>>(&found->second);
        if (radii == nullptr)
          throw ModelError("run record: telescope_r_m must be a float array");
        radii_ = *radii;
      }
    }
    event_.reset();
    writer_.write_metadata(obdeect::stream::RecordType::run, metadata);
    run_status_pending_ = !per_telescope();
  }

  // Sphere radius [m] of a telescope (id 1 is the first entry, as sim_telarray).
  [[nodiscard]] double sphere_radius(std::uint32_t telescope_id) const {
    if (telescope_id == 0 || telescope_id > radii_.size())
      throw ModelError("sphere run: telescope " + std::to_string(telescope_id) +
                       " has no radius in the run record (telescope_r_m)");
    return radii_[telescope_id - 1];
  }

  void efficiency_status(const Metadata &metadata) {
    using namespace stream_trace_detail;
    if (!telescope_) {
      for (auto &[name, factor] : factors_of(metadata))
        run_factors_.insert_or_assign(name, std::move(factor));
      // The run status of the shared optics is written again in the output (run_status_pending_):
      // an own factor that is already applied would be written twice, so it is refused here.
      if (shared_)
        for (const auto &name : shared_->own_factors()) {
          const auto found = run_factors_.find(name);
          if (found != run_factors_.end() && found->second.state == "applied")
            throw ModelError("the efficiency status says " + name +
                             " is already applied; trace-stream applies it in its trace and "
                             "cannot skip it");
        }
      bool normalise = shared_ && shared_->normalise();
      for (const auto &[id, optics] : per_telescope_)
        normalise = normalise || optics->normalise();
      if (normalise && contains(names_of(metadata), kConstantFactor))
        throw ModelError(std::string("the efficiency status already names ") + kConstantFactor);
      writer_.write_metadata(obdeect::stream::RecordType::efficiency_status, metadata);
      return;
    }
    Telescope &telescope = *telescope_;
    if (telescope.photons_seen)
      throw ModelError("telescope " + std::to_string(telescope.id) +
                       ": efficiency_status after its photons");
    if (!telescope.kept)
      return;
    for (auto &[name, factor] : factors_of(metadata))
      telescope.factors.insert_or_assign(name, std::move(factor));
    if (per_telescope()) {
      if (telescope.status_written)
        throw ModelError("telescope " + std::to_string(telescope.id) +
                         ": more than one efficiency_status before its photons");
      write_status(telescope, metadata);
      return;
    }
    writer_.write_metadata(obdeect::stream::RecordType::efficiency_status,
                           with_applied(metadata, optics(telescope.id)));
  }

  void event(const Metadata &metadata) {
    (void)stream_trace_detail::int_key(metadata, "event_number", "event");
    event_ = metadata;
    telescope_.reset();
    next_index_.clear();
    ++counts_.events;
    writer_.write_metadata(obdeect::stream::RecordType::event, metadata);
  }

  void telescope(const Metadata &metadata) {
    using namespace stream_trace_detail;
    const auto id = static_cast<std::uint32_t>(int_key(metadata, "telescope_id", "telescope"));
    Telescope telescope;
    telescope.id = id;
    telescope.kept = !settings_.telescopes || settings_.telescopes->count(id) != 0;
    const auto next = next_index_.find(id);
    telescope.out_index = next == next_index_.end() ? 0 : next->second;
    telescope_ = std::move(telescope);
    if (telescope_->kept && sphere_) {
      const double radius = sphere_radius(id);
      for (const auto &[name, distance] : optics(id).enclosure)
        if (distance > radius)
          throw ModelError("telescope " + std::to_string(id) +
                           ": the sphere (telescope_r_m = " + format_3f(radius) +
                           " m) must enclose the optical model, but '" + name + "' lies " +
                           format_3f(distance) + " m from the sphere centre; raise telescope_r_m");
    }
    if (telescope_->kept) {
      Metadata out = metadata;
      out["frame"] = static_cast<std::int64_t>(obdeect::stream::Frame::camera);
      writer_.write_metadata(obdeect::stream::RecordType::telescope, out);
    }
  }

  void write_status(Telescope &telescope, const Metadata &metadata) {
    writer_.write_metadata(obdeect::stream::RecordType::efficiency_status,
                           stream_trace_detail::own_status(optics(telescope.id), metadata));
    telescope.status_written = true;
  }

  void ensure_status(Telescope &telescope) {
    if (per_telescope() && !telescope.status_written)
      write_status(telescope, {});
  }

  void telescope_end(const Metadata &metadata) {
    using namespace stream_trace_detail;
    if (!telescope_ ||
        int_key(metadata, "telescope_id", "telescope_end") != std::int64_t{telescope_->id})
      throw ModelError("telescope_end record without its telescope record");
    Telescope &telescope = *telescope_;
    if (telescope.kept) {
      // Applies the efficiency status checks even for a telescope without photons.
      ensure_status(telescope);
      (void)config(telescope);
      Metadata out = metadata;
      out["photons"] = telescope.weight;
      writer_.write_metadata(obdeect::stream::RecordType::telescope_end, out);
      next_index_[telescope.id] = telescope.out_index;
    }
    telescope_.reset();
  }

  const stream_trace_detail::Config &config(Telescope &telescope) {
    using namespace stream_trace_detail;
    if (!telescope.config) {
      std::map<std::string, Factor> factors = run_factors_;
      for (const auto &[name, factor] : telescope.factors)
        factors.insert_or_assign(name, factor);
      Config config;
      for (const auto &name : optics(telescope.id).own_factors()) {
        const auto found = factors.find(name);
        if (found == factors.end() || found->second.state == "pending")
          continue;
        if (found->second.state == "applied")
          throw ModelError("telescope " + std::to_string(telescope.id) +
                           ": the efficiency status says " + name +
                           " is already applied; trace-stream applies it in its trace and "
                           "cannot skip it");
        config.envelopes.push_back(envelope_table(name, found->second));
      }
      config.roulette = config.envelopes.empty();
      telescope.config = std::move(config);
    }
    return *telescope.config;
  }

  // Pointing of the ground-frame photons: the settings, else the run record.
  [[nodiscard]] std::pair<double, double> pointing_for_run() const {
    if (settings_.pointing)
      return *settings_.pointing;
    if (run_pointing_)
      return *run_pointing_;
    throw ModelError("no telescope pointing for ground-frame photons: give "
                     "--pointing-zenith-deg and --pointing-azimuth-deg, or the run keys "
                     "telescope_pointing_zenith_rad and telescope_pointing_azimuth_rad");
  }

  // Input planes of a telescope's photons, cached per optical model and pointing (the same for
  // every event and block of the run).
  InputPlanes &input_planes(std::uint32_t telescope_id, bool ground_frame) {
    const Optics &model_optics = optics(telescope_id);
    std::optional<std::pair<double, double>> pointing;
    if (ground_frame)
      pointing = pointing_for_run();
    const GeometryKey key{&model_optics, pointing};
    auto found = geometry_.find(key);
    if (found == geometry_.end()) {
      found =
          geometry_
              .emplace(key, pointing ? stream_trace_detail::ground_input_planes(
                                           *model_optics.model, pointing->first, pointing->second)
                                     : stream_trace_detail::telescope_input_planes(model_optics))
              .first;
    }
    return found->second;
  }

  void photons(const obdeect::stream::PhotonChunk &chunk) {
    using namespace stream_trace_detail;
    using obdeect::stream::Column;
    if (!telescope_ || telescope_->id != chunk.telescope_id)
      throw ModelError("photons of telescope " + std::to_string(chunk.telescope_id) +
                       " outside their block");
    Telescope &telescope = *telescope_;
    if (!telescope.kept)
      return;
    if (!event_)
      throw ModelError("photons outside an event record");
    telescope.photons_seen = true;
    if (chunk.has(Column::p))
      throw ModelError("photon column p (acceptance already applied) is not supported");
    const std::size_t rows = chunk.rows();
    if (rows == 0)
      return;
    if (chunk.frame != obdeect::stream::Frame::ground &&
        chunk.frame != obdeect::stream::Frame::telescope)
      throw ModelError("photon frame " + std::to_string(static_cast<std::uint32_t>(chunk.frame)) +
                       " is not supported: trace-stream takes ground (1) or telescope (2) photons");
    // Sphere rows (G frame of a sphere run): each row starts at its own point, see below.
    const bool sphere_rows = sphere_ && chunk.frame == obdeect::stream::Frame::ground;
    if (sphere_rows) {
      try {
        st::validate_sphere_rows(chunk, sphere_radius(telescope.id));
      } catch (const std::invalid_argument &error) {
        throw ModelError(error.what());
      }
    }
    ensure_status(telescope);
    const Config &config = this->config(telescope);
    const Optics &optics = this->optics(telescope.id);
    const std::int64_t event_number = int_key(*event_, "event_number", "event");

    // Input bunches on the input plane. Photons moving down start on the upper plane, photons
    // moving up on the lower plane (below all optics); cz = 0 has no plane to reach.
    const auto &cx_in = chunk[Column::cx];
    const auto &cy_in = chunk[Column::cy];
    std::vector<double> cz(rows);
    bool cz_mismatch = false, bad_direction = false, horizontal = false, upward = false,
         downward = false;
    for (std::size_t i = 0; i < rows; ++i) {
      const double cx = cx_in[i], cy = cy_in[i];
      const double v = 1.0 - cx * cx - cy * cy;
      cz[i] = -std::sqrt(std::isnan(v) ? v : std::max(0.0, v));
      if (chunk.has(Column::cz)) {
        cz[i] = chunk[Column::cz][i];
        cz_mismatch = cz_mismatch || !st::cz_is_consistent(chunk[Column::cz][i], cx, cy);
      }
      bad_direction = bad_direction || !std::isfinite(cz[i]) || !st::direction_in_unit_disc(cx, cy);
      horizontal = horizontal || cz[i] == 0.0;
      upward = upward || cz[i] > 0.0;
      downward = downward || cz[i] < 0.0;
    }
    if (bad_direction)
      throw ModelError("photon directions must be finite unit vectors (cx, cy in the unit disc)");
    if (cz_mismatch)
      throw ModelError("photon column cz does not match cx and cy");
    // Plane rows cannot be horizontal (they have no plane to start on); sphere rows can.
    if (horizontal && !sphere_rows)
      throw ModelError("horizontal photons relative to the telescope axis are not supported yet");
    const bool ground = chunk.frame == obdeect::stream::Frame::ground;
    if (ground && chunk.has(Column::z) && !sphere_rows)
      for (const float z : chunk[Column::z])
        if (z != 0.0F)
          throw ModelError("ground-frame photons must lie in the plane z = 0");
    if (!ground && !chunk.has(Column::z))
      throw ModelError("telescope-frame photons need the column z");
    InputPlanes &input = input_planes(telescope.id, ground);
    const std::array<double, 9> rotation = input.rotation;
    const Vec offset = input.offset;
    // Sphere rows start at their own z: the input planes are not used.
    const double plane = sphere_rows ? 0.0 : input.upper;
    const double lower = sphere_rows ? std::numeric_limits<double>::quiet_NaN() : input.lower;

    if (event_number < 0 || event_number >= (std::int64_t{1} << 32))
      throw ModelError("event numbers must lie below 2^32 (the random key)");
    // Bunch ids are the photon indices of the block: the event is in the random key.
    if (chunk.first_index > st::max_photon_index ||
        rows - 1 > st::max_photon_index - chunk.first_index)
      throw ModelError("photon indices must lie below 2^40 (the random key)");
    if (downward && !sphere_rows && !input.upper_checked) {
      check_input_plane(*optics.model, rotation, offset, plane);
      input.upper_checked = true;
    }
    if (upward && !sphere_rows && !input.lower_checked) {
      check_lower_input_plane(*optics.model, rotation, offset, lower);
      input.lower_checked = true;
    }
    const auto &x_in = chunk[Column::x];
    const auto &y_in = chunk[Column::y];
    const auto &t_in = chunk[Column::t];
    std::vector<float> x(rows), y(rows), time(rows);
    std::vector<std::uint64_t> bunch_id(rows);
    // Path from the emission to each start point [m] (format 1.2 column path, carried through):
    // the input path, moved by the signed distance along the direction from the row to the start.
    const bool has_path = chunk.has(Column::path);
    std::vector<double> path_start(has_path ? rows : 0);
    for (std::size_t i = 0; i < rows; ++i) {
      if (sphere_rows) {
        // The row is the start point itself; its time and path are those there.
        x[i] = x_in[i];
        y[i] = y_in[i];
        time[i] = t_in[i];
      } else {
        const double z = ground ? 0.0 : static_cast<double>(chunk[Column::z][i]);
        const double row_plane = cz[i] < 0.0 ? plane : lower;
        const double s = (row_plane - z) / cz[i];
        x[i] = static_cast<float>(static_cast<double>(x_in[i]) + s * static_cast<double>(cx_in[i]));
        y[i] = static_cast<float>(static_cast<double>(y_in[i]) + s * static_cast<double>(cy_in[i]));
        const double norm = std::sqrt(static_cast<double>(cx_in[i]) * cx_in[i] +
                                      static_cast<double>(cy_in[i]) * cy_in[i] + cz[i] * cz[i]);
        time[i] = static_cast<float>(static_cast<double>(t_in[i]) + s * norm / kSpeedOfLightMPerNs);
        if (has_path)
          path_start[i] = static_cast<double>(chunk[Column::path][i]) + s * norm;
      }
      if (has_path && sphere_rows)
        path_start[i] = static_cast<double>(chunk[Column::path][i]);
      if (has_path && path_start[i] < 0.0)
        throw ModelError("photon path is shorter than its move to the start point");
      bunch_id[i] = chunk.first_index + i;
    }
    const std::vector<float> zeros(rows, 0.0F);
    const auto &zem = chunk.has(Column::zem) ? chunk[Column::zem] : zeros;
    const auto &wavelength = chunk.has(Column::wavelength) ? chunk[Column::wavelength] : zeros;
    const BunchBlock block{x.data(),
                           y.data(),
                           cx_in.data(),
                           cy_in.data(),
                           time.data(),
                           zem.data(),
                           chunk[Column::w].data(),
                           wavelength.data(),
                           bunch_id.data(),
                           rows,
                           chunk.has(Column::cz) ? chunk[Column::cz].data() : nullptr,
                           sphere_rows ? chunk[Column::z].data() : nullptr,
                           has_path ? path_start.data() : nullptr};
    BunchFrame frame;
    std::copy(rotation.begin(), rotation.end(), frame.rotation);
    std::copy(offset.begin(), offset.end(), frame.offset_m);
    frame.length_unit_m = 1.0;
    frame.plane_z = plane;
    frame.lower_plane_z = lower;
    BunchOptions options;
    options.roulette = config.roulette;
    options.normalise_transmission = optics.normalise();
    if (law_)
      options.wavelength_law = &*law_;
    std::vector<double> position, direction, t_out, weight, wavelength_out;
    std::vector<std::uint64_t> bunch_index;
    std::vector<std::uint32_t> sub_photon;
    std::vector<double> path_out;
    const std::size_t traced = trace_bunches(
        optics.model->model, block, frame, options,
        RandomKey{settings_.seed, static_cast<std::uint32_t>(event_number), telescope.id}, pool_,
        [&](std::size_t n) {
          position.resize(3 * n);
          direction.resize(3 * n);
          t_out.resize(n);
          weight.resize(n);
          wavelength_out.resize(n);
          bunch_index.resize(n);
          sub_photon.resize(n);
          if (has_path)
            path_out.resize(n);
          return BunchResult{position.data(),       direction.data(),
                             t_out.data(),          weight.data(),
                             wavelength_out.data(), bunch_index.data(),
                             sub_photon.data(),     has_path ? path_out.data() : nullptr};
        });
    counts_.photons_in += rows;

    // Camera-frame columns.
    if (!config.envelopes.empty()) {
      std::vector<double> inverse(traced, 1.0);
      for (const auto &table : config.envelopes) {
        if (traced > 0) {
          const auto [lo, hi] =
              std::minmax_element(wavelength_out.begin(), wavelength_out.begin() + traced);
          if (*lo < table.nodes.front() || *hi > table.nodes.back())
            throw ModelError("photon wavelengths lie outside the efficiency status envelope [" +
                             format_g(table.nodes.front()) + ", " + format_g(table.nodes.back()) +
                             "] nm");
        }
        for (std::size_t i = 0; i < traced; ++i)
          inverse[i] /= interpolate(wavelength_out[i], table.nodes, table.values);
      }
      for (std::size_t i = 0; i < traced; ++i)
        weight[i] = weight[i] * inverse[i];
    }
    const bool with_zem = chunk.has(Column::zem);
    for (std::size_t start = 0; start < traced; start += obdeect::stream::max_chunk_rows) {
      const std::size_t n = std::min<std::size_t>(obdeect::stream::max_chunk_rows, traced - start);
      std::uint32_t mask = 0;
      for (const Column c : {Column::x, Column::y, Column::cx, Column::cy, Column::cz, Column::t,
                             Column::w, Column::wavelength})
        mask |= obdeect::stream::column_bit(c);
      if (with_zem)
        mask |= obdeect::stream::column_bit(Column::zem);
      if (has_path)
        mask |= obdeect::stream::column_bit(Column::path);
      out_.resize(n, mask);
      out_.telescope_id = telescope.id;
      out_.frame = obdeect::stream::Frame::camera;
      out_.first_index = telescope.out_index;
      for (std::size_t k = 0; k < n; ++k) {
        const std::size_t i = start + k;
        out_[Column::x][k] = static_cast<float>(position[3 * i]);
        out_[Column::y][k] = static_cast<float>(position[3 * i + 1]);
        out_[Column::cx][k] = static_cast<float>(direction[3 * i]);
        out_[Column::cy][k] = static_cast<float>(direction[3 * i + 1]);
        out_[Column::cz][k] = static_cast<float>(direction[3 * i + 2]);
        out_[Column::t][k] = static_cast<float>(t_out[i]);
        out_[Column::w][k] = static_cast<float>(weight[i]);
        out_[Column::wavelength][k] = static_cast<float>(wavelength_out[i]);
        if (with_zem)
          out_[Column::zem][k] = zem[bunch_index[i]];
        if (has_path)
          out_[Column::path][k] = static_cast<float>(path_out[i]);
      }
      writer_.write_chunk(out_);
      telescope.out_index += n;
      // The sum of the weights in output order, so that it does not depend on the input chunking.
      for (const float w : out_[Column::w])
        telescope.weight += static_cast<double>(w);
    }
    counts_.photons_out += traced;
  }

  std::shared_ptr<const Optics> shared_;
  std::map<std::uint32_t, std::shared_ptr<const Optics>> per_telescope_;
  obdeect::stream::Writer &writer_;
  bool sphere_{false};        // The run's ground representation is sphere.
  std::vector<double> radii_; // telescope_r_m of a sphere run.
  StreamTraceSettings settings_;
  ThreadPool pool_;
  std::optional<std::pair<double, double>> run_pointing_;
  std::map<std::string, stream_trace_detail::Factor> run_factors_;
  std::optional<InverseSquareLaw> law_;
  std::optional<Metadata> event_;
  std::map<GeometryKey, InputPlanes> geometry_;
  std::optional<Telescope> telescope_;
  std::map<std::uint32_t, std::uint64_t> next_index_;
  StreamTraceCounts counts_;
  bool run_status_pending_{false};
  obdeect::stream::PhotonChunk out_;
};

// Traces a whole stream from `in` to `out`. `shared` is the model of every telescope, or null
// with `per_telescope` the model of each telescope id. Writes stream_end only after the whole
// input, so a failed run leaves a stream without it.
inline StreamTraceCounts trace_stream(
    std::istream &in, std::ostream &out, std::shared_ptr<const stream_trace_detail::Optics> shared,
    std::map<std::uint32_t, std::shared_ptr<const stream_trace_detail::Optics>> per_telescope,
    const StreamTraceSettings &settings) {
  namespace st = obdeect::stream;
  st::Reader reader(in);
  if (settings.telescopes) // Photon chunks of other telescopes are skipped, not decoded.
    reader.keep_photons_of({settings.telescopes->begin(), settings.telescopes->end()});
  st::Record record;
  if (!reader.next(record) || record.type != st::RecordType::stream)
    throw ModelError("input is not a photon stream");
  std::string producer;
  if (const auto found = record.metadata.find("producer"); found != record.metadata.end())
    producer = stream_trace_detail::value_text(found->second);
  st::Writer writer(out, producer.empty() ? "obdeect " + settings.version
                                          : producer + "; obdeect " + settings.version);
  StreamTracer tracer(std::move(shared), std::move(per_telescope), writer, settings);
  while (reader.next(record))
    tracer.handle(record);
  if (!reader.ended())
    throw std::runtime_error("truncated photon stream (no stream_end record)");
  if (!reader.at_clean_end())
    throw std::runtime_error("data after stream_end");
  writer.finish();
  return tracer.counts();
}

} // namespace obdeect::ray_tracing
