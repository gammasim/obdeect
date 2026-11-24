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
  const auto actual = trace_photon(with_window, ray, 400.0, 0.0, 1.0, 0, 1,
                                   {.apply_responses = true, .direct = true});
  const auto reference = trace_photon(without_window, ray, 400.0, 0.0, 1.0, 0, 1,
                                      {.apply_responses = true, .direct = true});
  CHECK(actual.status == Status::detected);
  const double perpendicular_shift =
      std::abs(actual.position.x - reference.position.x) * std::cos(theta);
  const double expected_shift =
      thickness * std::sin(theta) *
      (1.0 - std::cos(theta) / std::sqrt(n * n - std::sin(theta) * std::sin(theta)));
  CHECK_NEAR(perpendicular_shift, expected_shift, 1e-12);

  const auto normal = trace_photon(with_window, {{0, 0, -1}, {0, 0, 1}}, 400, 0, 1, 0, 1,
                                   {.apply_responses = true, .direct = true});
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
                                   0, 1, 0, 1, {.apply_responses = true, .direct = true});
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
                                 1, 0, 1, {.apply_responses = true, .direct = true});
  CHECK(lost.status == Status::camera_window_total_reflection);
}

void rim_and_concentric_sphere() {
  auto flat = direct_model();
  const auto rim = trace_photon(flat, {{1.5, 0, -1}, {0, 0, 1}}, 400, 0, 1, 0, 1,
                                {.apply_responses = true, .direct = true});
  CHECK(rim.status == Status::blocked_camera_window);
  CHECK(rim.terminal_surface == flat.camera_window->id);
  // Enters inside the aperture (r = 0.99 m), leaves outside it (r = 1.01 m).
  const auto exit_rim =
      trace_photon(flat, {{0.99 - std::tan(0.3), 0, -1}, {std::sin(0.3), 0, std::cos(0.3)}}, 400, 0,
                   1, 0, 1, {.apply_responses = true, .direct = true});
  CHECK(exit_rim.status == Status::blocked_camera_window);
  CHECK(exit_rim.position.x > 1.0);

  auto sphere = direct_model(-0.2, 0.1, 1.5, 1.0);
  sphere.camera_window->shape = Window::Shape::sphere;
  sphere.camera_window->radius_of_curvature = 3.5;
  const auto radial = trace_photon(sphere, {{0, 0, -1}, {0, 0, 1}}, 400, 0, 1, 0, 1,
                                   {.apply_responses = true, .direct = true});
  CHECK(radial.status == Status::detected);
  CHECK_NEAR(radial.position.x, 0.0, 1e-12);
  CHECK_NEAR(radial.position.y, 0.0, 1e-12);
  CHECK_NEAR(radial.direction.x, 0.0, 1e-12);
  CHECK_NEAR(radial.direction.z, 1.0, 1e-12);
}

// LST-like camera window: spherical outer surface (R = 3.5 m, 5 mm glass, n = 1.5, tau = 1),
// focal plane at z = 28 m, 20 mm flange between the window apex and the focal plane.
OpticalModel lst_like_model() {
  constexpr double radius = 3.5;
  constexpr double aperture = 1.17575;
  constexpr double focal_z = 28.0;
  constexpr double flange = 0.02;
  const double sag = radius - std::sqrt(radius * radius - aperture * aperture);
  auto model = direct_model(focal_z - flange - sag, 0.005, 1.5, 1.0);
  model.focal_surface = {1, focal_z, std::nullopt};
  Window &window = *model.camera_window;
  window.shape = Window::Shape::sphere;
  window.radius_of_curvature = radius;
  window.aperture_radius = aperture;
  window.outer_radius = 1.2517;
  return model;
}

// A photon between the window and the focal plane is inside the outer sphere. Its first
// window root is behind it (the far side of the sphere), so it must not be refracted there.
void photon_behind_window_is_not_refracted() {
  const auto model = lst_like_model();
  const auto behind = trace_photon(model, {{0.1, 0.0, 27.95}, {0.0, 0.0, 1.0}}, 400.0, 0.0, 1.0, 0,
                                   1, {.apply_responses = true, .direct = true});
  CHECK(behind.status == Status::detected);
  CHECK_NEAR(behind.position.z, 28.0, 1e-12);
  CHECK_NEAR(behind.position.x, 0.1, 1e-12);
  CHECK_NEAR(behind.weight, 1.0, 1e-12);

  // The same direction from in front of the window does cross it (refraction lowers the weight).
  const auto front = trace_photon(model, {{0.1, 0.0, 27.0}, {0.0, 0.0, 1.0}}, 400.0, 0.0, 1.0, 0, 1,
                                  {.apply_responses = true, .direct = true});
  CHECK(front.status == Status::detected);
  CHECK(front.weight < 1.0);
}

// Flat plate: a parallel plate keeps the direction and shifts the axial crossing of the
// focus by dz = t * (1 - cos(theta) / sqrt(n^2 - sin^2(theta))).
void flat_plate_focus_shift() {
  constexpr double n = 1.5;
  constexpr double thickness = 0.005;
  auto model = direct_model(0.0, thickness, n, 1.0);
  model.camera_window->aperture_radius = 10.0;
  model.camera_window->outer_radius = 11.0;
  const auto shift_at = [&model](double theta) {
    const double s = std::sin(theta);
    const double c = std::cos(theta);
    const Ray ray{{-2.0 * std::tan(theta), 0.0, -1.0}, {s, 0.0, c}};
    const auto result =
        trace_photon(model, ray, 400.0, 0.0, 1.0, 0, 1, {.apply_responses = true, .direct = true});
    CHECK(result.status == Status::detected);
    // Without the plate the ray would reach x = 0 on the focal plane (z = 1 m).
    return -result.position.x / std::tan(theta);
  };
  for (const double theta : {0.05, 0.2, 0.35}) {
    const double s = std::sin(theta);
    const double c = std::cos(theta);
    CHECK_NEAR(shift_at(theta), thickness * (1.0 - c / std::sqrt(n * n - s * s)), 1e-12);
  }
  CHECK_NEAR(shift_at(1e-4), thickness * (1.0 - 1.0 / n), 1e-9);
}

// Rays aimed at the focal point through the LST-like shell. Reference arrival x on the focal
// plane (m) from an independent numpy vector ray trace (sphere intersections and vector Snell,
// written separately from this code, 2026-11-24).
void curved_shell_focus_shift_lst_geometry() {
  const auto model = lst_like_model();
  constexpr double reference[][2] = {
      {0.05, -0.07327860303215673e-3},
      {0.2, -0.30159248778743875e-3},
      {0.35, -0.5619972912839222e-3},
  };
  for (const auto &[theta, expected_x] : reference) {
    const Ray ray{{-std::tan(theta), 0.0, 27.0}, {std::sin(theta), 0.0, std::cos(theta)}};
    const auto result =
        trace_photon(model, ray, 400.0, 0.0, 1.0, 0, 1, {.apply_responses = true, .direct = true});
    CHECK(result.status == Status::detected);
    CHECK_NEAR(result.position.z, 28.0, 1e-12);
    CHECK_NEAR(result.position.x, expected_x, 1e-12);
  }

  // A ray aimed at the common centre of curvature meets both surfaces at normal incidence.
  const double sag = 3.5 - std::sqrt(3.5 * 3.5 - 1.17575 * 1.17575);
  const double centre_z = 28.0 - 0.02 - sag + 3.5;
  const Vec3 origin{0.5, 0.0, 27.0};
  const Vec3 toward{0.0 - origin.x, 0.0, centre_z - origin.z};
  const double length = std::sqrt(dot(toward, toward));
  const Vec3 direction{toward.x / length, 0.0, toward.z / length};
  const auto radial = trace_photon(model, {origin, direction}, 400.0, 0.0, 1.0, 0, 1,
                                   {.apply_responses = true, .direct = true});
  const double n = 1.5;
  const double t0 = 1.0 - std::pow((n - 1.0) / (n + 1.0), 2);
  CHECK(radial.status == Status::detected);
  CHECK_NEAR(radial.direction.x, direction.x, 1e-12);
  CHECK_NEAR(radial.direction.z, direction.z, 1e-12);
  CHECK_NEAR(radial.weight, t0 * t0, 1e-12);
}

// The document with its sha256 entry set to the hash of its content (load_model checks it).
std::string signed_document(std::string text) {
  const std::string placeholder = "\"sha256\":\"ignored\"";
  const auto at = text.find(placeholder);
  if (at != std::string::npos)
    text.replace(at, placeholder.size(), "\"sha256\":\"" + model_hash(text) + "\"");
  return text;
}

void json_window_round_trip() {
  const std::string doc =
      R"({"format":"obdeect.optical-model.v2","sha256":"ignored","sampling":{"disk_radius_m":1.2,"disk_z_m":0,"start_z_m":2},
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
  const auto model = load_model(signed_document(doc));
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
  CHECK(load_model(signed_document(without_window)).camera_window == std::nullopt);
  CHECK(model_hash(doc) != model_hash(without_window));
}

} // namespace

int main() {
  // run_checks reports CHECK failures in the exit status (and catches an escaping exception).
  return run_checks([] {
    flat_plate_shift_and_weight();
    fresnel_at_sixty_degrees_and_total_reflection();
    rim_and_concentric_sphere();
    photon_behind_window_is_not_refracted();
    flat_plate_focus_shift();
    curved_shell_focus_shift_lst_geometry();
    json_window_round_trip();
  });
}
