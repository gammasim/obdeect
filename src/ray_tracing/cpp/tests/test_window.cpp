// Camera-window geometry, Fresnel transport, and optical-model loading.
#include "check.hpp"

#include "obdeect/ray_tracing/model_hash.hpp"
#include "obdeect/ray_tracing/model_json.hpp"
#include "obdeect/ray_tracing/trace.hpp"

#include <cmath>
#include <string>

using namespace obdeect::ray_tracing;

namespace {

Response constant(double value) { return {300.0, 1.0, 1, {0.0}, {value}}; }

OpticalModel direct_model(double z = 0.0, double thickness = 0.1, double n = 1.5,
                          double tau = 0.9) {
  OpticalModel model;
  model.focal_surface = {1, 1.0, std::nullopt};
  Window window;
  window.id = 2;
  window.outer_apex_z = z;
  window.thickness = thickness;
  window.aperture_radius = 1.0;
  window.outer_radius = 2.0;
  window.refractive_index = constant(n);
  window.internal_transmittance = constant(tau);
  model.camera_window = window;
  return model;
}

double fresnel_transmission(double n1, double n2, double cos_i) {
  const double sin2_t = (n1 / n2) * (n1 / n2) * (1.0 - cos_i * cos_i);
  const double cos_t = std::sqrt(1.0 - sin2_t);
  const double rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t);
  const double rp = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t);
  return 1.0 - 0.5 * (rs * rs + rp * rp);
}

void flat_plate_shift_and_weight() {
  constexpr double theta = 0.3;
  constexpr double n = 1.5;
  constexpr double thickness = 0.1;
  constexpr double tau = 0.9;
  auto with_window = direct_model(0.0, thickness, n, tau);
  auto without_window = with_window;
  without_window.camera_window.reset();
  const Ray ray{{0.0, 0.0, -1.0}, {std::sin(theta), 0.0, std::cos(theta)}};
  const auto actual = trace_photon(with_window, ray, 400.0, 0.0, 1.0, 0, 1, true, true);
  const auto reference = trace_photon(without_window, ray, 400.0, 0.0, 1.0, 0, 1, true, true);
  CHECK(actual.status == Status::detected);
  const double perpendicular_shift =
      std::abs(actual.position.x - reference.position.x) * std::cos(theta);
  const double expected_shift =
      thickness * std::sin(theta) *
      (1.0 - std::cos(theta) / std::sqrt(n * n - std::sin(theta) * std::sin(theta)));
  CHECK_NEAR(perpendicular_shift, expected_shift, 1e-12);

  const auto normal =
      trace_photon(with_window, {{0, 0, -1}, {0, 0, 1}}, 400, 0, 1, 0, 1, true, true);
  const double t0 = (1.0 - std::pow((n - 1.0) / (n + 1.0), 2));
  CHECK_NEAR(normal.weight, t0 * t0 * tau, 1e-12);
  CHECK_NEAR(normal.time_ns - 2.0 / kSpeedOfLightMPerNs,
             (n - 1.0) * thickness / kSpeedOfLightMPerNs, 1e-12);
}

void fresnel_at_sixty_degrees_and_total_reflection() {
  constexpr double theta = 1.04719755119659774615;
  constexpr double n = 1.5;
  constexpr double tau = 0.9;
  auto model = direct_model(0, 0.1, n, tau);
  model.camera_window->aperture_radius = 10.0;
  model.camera_window->outer_radius = 11.0;
  const auto result = trace_photon(model, {{0, 0, -1}, {std::sin(theta), 0, std::cos(theta)}}, 400,
                                   0, 1, 0, 1, true, true);
  const double sin_t = std::sin(theta) / n;
  const double cos_t = std::sqrt(1.0 - sin_t * sin_t);
  const double t0 = fresnel_transmission(1.0, n, std::cos(theta));
  const double t1 = fresnel_transmission(n, 1.0, cos_t);
  CHECK(result.status == Status::detected);
  CHECK_NEAR(result.weight, t0 * t1 * std::pow(tau, 1.0 / cos_t), 1e-12);

  auto tir = direct_model(0, 0.1, 0.5, tau);
  tir.camera_window->aperture_radius = 10.0;
  tir.camera_window->outer_radius = 11.0;
  const auto lost = trace_photon(tir, {{0, 0, -1}, {std::sin(theta), 0, std::cos(theta)}}, 400, 0,
                                 1, 0, 1, true, true);
  CHECK(lost.status == Status::camera_window_total_reflection);
}

void rim_and_concentric_sphere() {
  auto flat = direct_model();
  const auto rim = trace_photon(flat, {{1.5, 0, -1}, {0, 0, 1}}, 400, 0, 1, 0, 1, true, true);
  CHECK(rim.status == Status::blocked_camera_window);
  CHECK(rim.terminal_surface == flat.camera_window->id);
  // Enters inside the aperture (r = 0.99 m), leaves outside it (r = 1.01 m).
  const auto exit_rim =
      trace_photon(flat, {{0.99 - std::tan(0.3), 0, -1}, {std::sin(0.3), 0, std::cos(0.3)}}, 400, 0,
                   1, 0, 1, true, true);
  CHECK(exit_rim.status == Status::blocked_camera_window);
  CHECK(exit_rim.position.x > 1.0);

  auto sphere = direct_model(-0.2, 0.1, 1.5, 1.0);
  sphere.camera_window->shape = Window::Shape::sphere;
  sphere.camera_window->radius_of_curvature = 3.5;
  const auto radial = trace_photon(sphere, {{0, 0, -1}, {0, 0, 1}}, 400, 0, 1, 0, 1, true, true);
  CHECK(radial.status == Status::detected);
  CHECK_NEAR(radial.position.x, 0.0, 1e-12);
  CHECK_NEAR(radial.position.y, 0.0, 1e-12);
  CHECK_NEAR(radial.direction.x, 0.0, 1e-12);
  CHECK_NEAR(radial.direction.z, 1.0, 1e-12);
}

void json_window_round_trip() {
  const std::string doc = R"({"format":"obdeect.optical-model.v2","sha256":"ignored",
    "provenance":{},"notes":[],"optics":{"propagation_group_index":1.0,
    "primary":{"facets":[{"id":0,"centre_m":[0,0,0],"normal":[0,0,1],
    "tangent":[1,0,0],"shape":"circle","diameter_m":1,"radius_of_curvature_m":0}]},
    "secondary":null,"secondary_blocks_incoming":false,"focal_surface":{"id":1,"plane_z_m":1},
    "camera_window":{"id":2,"shape":"sphere","outer_apex_z_m":0,
    "radius_of_curvature_m":3.5,"thickness_m":0.005,"aperture_radius_m":1.17575,
    "outer_radius_m":1.2517,"refractive_index":{"wavelength_min_nm":300,
    "wavelength_step_nm":1,"incidence_deg":[0],"values":[[1.49],[1.5]]},
    "internal_transmittance":{"wavelength_min_nm":300,"wavelength_step_nm":1,
    "incidence_deg":[0],"values":[[0.9],[0.8]]}},"obscurers":[]}})";
  const auto model = load_model(doc);
  CHECK(model.camera_window.has_value());
  CHECK(model.camera_window->shape == Window::Shape::sphere);
  CHECK_NEAR(model.camera_window->radius_of_curvature, 3.5, 1e-15);
  auto without_window = doc;
  const auto first = without_window.find("\"camera_window\":{");
  const auto preceding_comma = without_window.rfind(',', first);
  const auto following_comma = without_window.find(",\"obscurers\":[]", first);
  CHECK(first != std::string::npos && preceding_comma != std::string::npos &&
        following_comma != std::string::npos);
  without_window.erase(preceding_comma, following_comma - preceding_comma);
  CHECK(load_model(without_window).camera_window == std::nullopt);
  CHECK(model_hash(doc) != model_hash(without_window));
}

} // namespace

int main() {
  flat_plate_shift_and_weight();
  fresnel_at_sixty_degrees_and_total_reflection();
  rim_and_concentric_sphere();
  json_window_round_trip();
}
