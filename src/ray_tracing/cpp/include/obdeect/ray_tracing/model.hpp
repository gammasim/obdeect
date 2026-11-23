#pragma once

#include "obdeect/ray_tracing/obscurers.hpp"
#include "obdeect/ray_tracing/surfaces.hpp"
#include "obdeect/ray_tracing/table.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace obdeect::ray_tracing {

// Thrown for any invalid or unsupported optical model; for a model file the
// message names the JSON path.
struct ModelError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Random reflection error: with probability fraction2 the Gaussian width is
// sigma2, otherwise sigma1 (radians, per axis). `surface_slopes` perturbs the
// surface slope by half the angle; otherwise the reflected direction is rotated.
struct Scatter {
  double sigma1_rad{}, fraction2{}, sigma2_rad{};
  bool surface_slopes{};
};

struct Mirror {
  std::vector<Facet> facets; // Segmented mirror, or ...
  FacetGrid grid;
  std::optional<Asphere> asphere; // ... a continuous asphere with optional segments.
  std::vector<Segment> segments;
  std::vector<Segment::TangentFrame> segment_frames; // One per segment, from the asphere.
  SegmentGrid segment_grid;                          // Segments by x-y cell, sorted by id.
  std::optional<Response> reflectivity;
  std::optional<Scatter> scatter;
  // Relative reflectivity map, in facet map coordinates (Facet::map_u,
  // map_v) for facets and in telescope x-y for an asphere.
  std::optional<Map2D> degradation;
  // Facet mirrors only: z of the dish vertex in the telescope frame [m] (the optical axis point of
  // the dish sphere, at -mirror_offset). Optional; used for the in-telescope path. Not used in
  // tracing.
  std::optional<double> vertex_z_m;
};

// Direction-dependent transmission T = on_axis / (1 + amplitude (sin(theta) /
// scale)^power)^outer_power, theta measured from the optical axis.
struct TelescopeTransmission {
  double on_axis{1.0}, amplitude{}, angular_scale_rad{1.0}, power{2.0}, outer_power{1.0};

  [[nodiscard]] double at(const Vec3 &direction) const {
    // Constant transmission: exact, as pow(1, p) == 1 and on_axis / 1.0 == on_axis.
    if (amplitude == 0.0)
      return on_axis;
    const double s = std::hypot(direction.x, direction.y);
    return on_axis /
           std::pow(1.0 + amplitude * std::pow(s / angular_scale_rad, power), outer_power);
  }
};

// Photons end on the focal surface: a plane z = plane_z or an asphere. The
// plane is unbounded, so the image is never truncated.
struct FocalSurface {
  std::uint32_t id{};
  double plane_z{};
  std::optional<Asphere> asphere;
};

// Optional refractive camera window. Spherical radii are measured from the
// common centre; the outer sphere's vertex is at outer_apex_z.
struct Window {
  enum class Shape : std::uint8_t { flat, sphere };
  std::uint32_t id{};
  Shape shape{Shape::flat};
  double outer_apex_z{};
  double radius_of_curvature{};
  double thickness{};
  double aperture_radius{};
  double outer_radius{};
  Response refractive_index;
  Response internal_transmittance;
};

struct OpticalModel {
  Mirror primary;
  std::optional<Mirror> secondary;
  // Without a dedicated shadow, light from the sky that meets the back of the
  // secondary mirror is absorbed there.
  bool secondary_blocks_incoming{};
  FocalSurface focal_surface;
  std::optional<Window> camera_window;
  Obscurers obscurers;
  std::optional<TelescopeTransmission> telescope_transmission;
  double propagation_group_index{1.0};
  // Optical origin offset O [m, telescope frame]: a ground point P is traced at R P - O.
  std::array<double, 3> axes_offset_m{};
};

} // namespace obdeect::ray_tracing
