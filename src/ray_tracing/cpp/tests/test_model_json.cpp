// Optical-model loading: valid documents load; every error names its path.
#include "check.hpp"

#include "obdeect/ray_tracing/bunches.hpp"
#include "obdeect/ray_tracing/model_json.hpp"

#include <algorithm>
#include <array>
#include <clocale>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace obdeect::ray_tracing;

namespace {

const std::string kValid = R"({
  "format": "obdeect.optical-model.v2",
  "sha256": "ignored by the loader",
  "provenance": {"model": "test"},
  "notes": [],
  "optics": {
    "propagation_group_index": 1.0003,
    "telescope_transmission": {"on_axis": 0.9, "amplitude": 0.0, "angular_scale_rad": 0.1,
                               "power": 2.0, "outer_power": 1.0},
    "primary": {
      "facets": [{"id": 0, "centre_m": [0, 0, 0], "normal": [0, 0, 1], "tangent": [1, 0, 0],
                  "shape": "hexagon_flat_y", "diameter_m": 1.2, "radius_of_curvature_m": 32.0}],
      "reflectivity": {"wavelength_min_nm": 300, "wavelength_step_nm": 1, "incidence_deg": [0],
                       "values": [[0.8], [0.9]]},
      "scatter": {"sigma1_rad": 0.0001, "fraction2": 0.2, "sigma2_rad": 0.0004,
                  "method": "outgoing_angles"},
      "degradation": {"x_m": [-1, 1], "y_m": [-1, 1], "values": [[1, 1], [1, 0.9]], "outside": "clamp"}
    },
    "secondary": null,
    "secondary_blocks_incoming": false,
    "focal_surface": {"id": 1, "plane_z_m": 16.0},
    "obscurers": [
      {"id": 2, "name": "mast-0", "group": "camera_support", "legs": ["incoming", "after_primary"],
       "shape": "cylinder", "start_m": [3, 0, 0], "end_m": [0, 0, 16], "diameter_m": 0.2},
      {"id": 3, "name": "box", "group": "dish_structure", "legs": ["incoming"], "shape": "polyhedron",
       "vertices_m": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
       "planes": [[0, 0, -1, 0], [0, -1, 0, 0], [-1, 0, 0, 0],
                  [0.5773502691896258, 0.5773502691896258, 0.5773502691896258, 0.5773502691896258]]},
      {"id": 4, "name": "camera", "group": "camera", "legs": ["incoming"], "shape": "plate",
       "centre_m": [0, 0, 16], "normal": [0, 0, 1], "tangent": [1, 0, 0], "aperture": "square",
       "diameter_m": 2.0}
    ]
  }
})";

// Load `text` with `from` replaced by `to`; return the error message, if any.
std::string error_after(const std::string &from, const std::string &to) {
  std::string text = kValid;
  const auto at = text.find(from);
  if (at == std::string::npos)
    return "test error: '" + from + "' not found";
  text.replace(at, from.size(), to);
  try {
    (void)load_model(text);
  } catch (const ModelError &error) {
    return error.what();
  }
  return "";
}

bool contains(const std::string &text, const std::string &part) {
  return text.find(part) != std::string::npos;
}

// Columns of one trace_bunches call without roulette: positions and their row count.
struct Traced {
  std::vector<double> position;
  std::size_t rows{};
};

Traced trace_bunches_of(const OpticalModel &model, const BunchBlock &in, const BunchFrame &frame) {
  Traced out;
  ThreadPool pool(2);
  std::vector<double> direction, time, weight, wavelength;
  std::vector<std::uint64_t> bunch;
  std::vector<std::uint32_t> sub;
  const auto allocate = [&](std::size_t rows) {
    out.position.resize(3 * rows);
    direction.resize(3 * rows);
    time.resize(rows);
    weight.resize(rows);
    wavelength.resize(rows);
    bunch.resize(rows);
    sub.resize(rows);
    return BunchResult{out.position.data(), direction.data(), time.data(), weight.data(),
                       wavelength.data(),   bunch.data(),     sub.data()};
  };
  out.rows = trace_bunches(model, in, frame, BunchOptions{}, 7, pool, allocate);
  return out;
}

// The loaded axes offset is applied by the tracer: a ground point P is traced at R P - O, the
// same as the model without offset and frame offset -O.
void axes_offset_is_applied_by_the_tracer() {
  const std::string key = "\"optics\": {";
  std::string text = kValid;
  text.insert(text.find(key) + key.size(), "\n    \"axes_offset_m\": [0.5, -0.2, 0.1],");
  auto shifted = load_model(text);
  auto plain = load_model(kValid);
  CHECK((shifted.axes_offset_m == std::array<double, 3>{0.5, -0.2, 0.1}));
  CHECK((plain.axes_offset_m == std::array<double, 3>{}));
  // Without the obscurers, which would shadow the photons, the photons reach the primary.
  shifted.obscurers = Obscurers(std::vector<Obscurer>{}, std::vector<Plane>{});
  plain.obscurers = Obscurers(std::vector<Obscurer>{}, std::vector<Plane>{});

  constexpr std::size_t n = 200;
  std::vector<float> x(n), y(n), zero(n, 0.0F), size(n, 1.0F), wavelength(n, 400.0F), time(n, 0.0F);
  std::vector<std::uint64_t> ids(n);
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = static_cast<float>(0.2 + 0.6 * static_cast<double>(i) / n);
    y[i] = static_cast<float>(-0.4 + 0.8 * static_cast<double>((37 * i) % n) / n);
    ids[i] = i;
  }
  const BunchBlock in{x.data(),    y.data(),    zero.data(),       zero.data(), time.data(),
                      zero.data(), size.data(), wavelength.data(), ids.data(),  n};
  BunchFrame origin;
  origin.plane_z = 30.0;
  BunchFrame explicit_offset = origin;
  explicit_offset.offset_m[0] = -0.5;
  explicit_offset.offset_m[1] = 0.2;
  explicit_offset.offset_m[2] = -0.1;

  const Traced offset_model = trace_bunches_of(shifted, in, origin);
  const Traced reference = trace_bunches_of(plain, in, explicit_offset);
  CHECK(offset_model.rows > 0);
  CHECK(offset_model.rows == reference.rows);
  CHECK(offset_model.position.size() == reference.position.size());
  double largest = 0.0;
  for (std::size_t i = 0; i < std::min(offset_model.position.size(), reference.position.size());
       ++i)
    largest = std::max(largest, std::abs(offset_model.position[i] - reference.position[i]));
  CHECK(largest <= 1e-12);
  // The offset is not a no-op: the same model without it traces the photons elsewhere.
  const Traced unshifted = trace_bunches_of(plain, in, origin);
  CHECK(unshifted.position != offset_model.position);
}

// Numbers parse the same under a locale with a decimal comma. Without one installed the
// test cannot exercise the path: it is skipped, unless OBDEECT_RAY_TRACING_REQUIRE_COMMA_LOCALE is
// set (CI does), in which case a missing locale fails.
void decimal_comma_locale() {
  const std::string previous = std::setlocale(LC_NUMERIC, nullptr);
  bool comma = false;
  for (const char *name : {"de_DE.UTF-8", "de_DE.utf8", "de_DE", "fr_FR.UTF-8", "fr_FR.utf8"})
    if (std::setlocale(LC_NUMERIC, name)) {
      comma = true;
      break;
    }
  if (!comma) {
    CHECK(std::getenv("OBDEECT_RAY_TRACING_REQUIRE_COMMA_LOCALE") == nullptr);
    std::puts("decimal_comma_locale: no decimal-comma locale installed, skipped");
    return;
  }
  const auto value = json::Parser::parse(R"([1.5, -2.25e-3, 12345.678901234567, 0.1])");
  std::setlocale(LC_NUMERIC, previous.c_str());
  const auto &numbers = std::get<json::Value::Array>(value.data);
  CHECK(numbers.size() == 4);
  CHECK(std::get<double>(numbers[0].data) == 1.5);
  CHECK(std::get<double>(numbers[1].data) == -2.25e-3);
  CHECK(std::get<double>(numbers[2].data) == 12345.678901234567);
  CHECK(std::get<double>(numbers[3].data) == 0.1);
}

} // namespace

void all_tests() {
  decimal_comma_locale();
  const auto model = load_model(kValid);
  CHECK(model.primary.facets.size() == 1);
  CHECK_NEAR(model.primary.facets[0].half_width, 0.6, 0.0);
  CHECK(model.primary.reflectivity && model.primary.reflectivity->wavelength_count == 2);
  CHECK(model.obscurers.items().size() == 3 && model.obscurers.planes().size() == 4);
  CHECK(model.obscurers.items()[0].legs == (Leg::incoming | Leg::after_primary));
  CHECK(!model.secondary && model.telescope_transmission);

  CHECK(contains(
      error_after(R"("normal": [0, 0, 1], "tangent")", R"("normal": [0, 0, 2], "tangent")"),
      "model.optics.primary.facets[0].normal: expected a unit vector"));
  CHECK(contains(error_after(R"("format": "obdeect.optical-model.v2")", R"("format": "v1")"),
                 "model.format"));
  CHECK(contains(error_after(R"("notes": [],)", R"("notes": [], "extra": 1,)"),
                 "unknown key 'extra'"));
  CHECK(contains(error_after(R"("id": 3, "name": "box")", R"("id": 2, "name": "box")"),
                 "used twice"));
  CHECK(contains(error_after(R"("values": [[0.8], [0.9]])", R"("values": [[0.8], [1.9]])"),
                 "reflectivity.values[1][0]: expected a value in [0, 1]"));
  CHECK(contains(error_after(R"("legs": ["incoming"], "shape": "plate")",
                             R"("legs": ["sideways"], "shape": "plate")"),
                 "obscurers[2].legs[0]: unknown leg"));
  CHECK(contains(error_after("0.5773502691896258]]", "0.1]]"), "vertex lies outside a face plane"));
  CHECK(contains(error_after(R"("plane_z_m": 16.0)", R"("asphere_x": 1)"),
                 "focal_surface: unknown key"));
  CHECK(contains(error_after("\"optics\": {", "\"optics\": {,"), "JSON syntax error"));
  CHECK(contains(error_after("\"optics\": {", "\"optics\": {\"axes_offset_m\": [1, 2],"),
                 "expected three components"));
  axes_offset_is_applied_by_the_tracer();
}

int main() { return run_checks(all_tests); }
