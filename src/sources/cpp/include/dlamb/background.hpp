#pragma once

// Night-sky background: photons start on the sphere of one telescope and travel inwards, with
// directions uniform in solid angle inside a cone (docs/spec/background.md).

#include "dlamb/density.hpp"
#include "dlamb/source.hpp"
#include "dlamb/vec3.hpp"

#include "gstrickd/photons.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace dlamb {

struct BackgroundSource {
  Vec3 centre{};
  double radius_m{};
  Vec3 pointing{0.0, 0.0, -1.0};
  double viewcone_rad{pi / 2.0};
  Spectrum spectrum{Fixed{400.0}};
  double pulse_width_ns{};
  double pulse_delay_ns{};
};

// Rotation about the y axis and about the z axis.
inline Vec3 rotate_y(const Vec3 &v, double a) {
  return {v.x * std::cos(a) + v.z * std::sin(a), v.y, v.z * std::cos(a) - v.x * std::sin(a)};
}
inline Vec3 rotate_z(const Vec3 &v, double a) {
  return {v.x * std::cos(a) - v.y * std::sin(a), v.x * std::sin(a) + v.y * std::cos(a), v.z};
}

inline std::uint64_t background_count(double mean_bunches, std::uint64_t seed, std::uint32_t event,
                                      std::uint32_t telescope) {
  return round_count(mean_bunches, seed, event, telescope);
}

inline Photon sample_background(const BackgroundSource &s, std::uint64_t seed, std::uint32_t event,
                                std::uint32_t telescope, std::uint64_t index) {
  const Uniforms ud = philox_draw(seed, event, telescope, index, Stream::direction);
  const Vec3 axis = normalised(s.pointing);
  const Basis basis = orthonormal_basis(axis);
  const double cos_theta = 1.0 - ud[0] * (1.0 - std::cos(s.viewcone_rad));
  const double sin_theta = std::sqrt(std::max(0.0, (1.0 - cos_theta) * (1.0 + cos_theta)));
  const Vec3 d = direction_about(axis, basis, cos_theta, sin_theta, two_pi * ud[1]);

  const Uniforms ut = philox_draw(seed, event, telescope, index, Stream::target);
  const double ri = 1.0 - std::sqrt(ut[0]);
  const double phi = two_pi * ut[1];
  const double zi = std::sqrt(std::max(0.0, 1.0 - ri * ri));
  const Vec3 local{ri * std::cos(phi), ri * std::sin(phi), zi};
  const double za = std::acos(std::clamp(-d.z, -1.0, 1.0));
  const double az = std::atan2(-d.y, -d.x);
  const Vec3 unit = rotate_z(rotate_y(local, za), az);

  Photon p;
  p.position = s.centre + unit * s.radius_m;
  p.direction = d;
  double jitter = 0.0;
  if (s.pulse_width_ns > 0.0) {
    const double u = philox_draw(seed, event, telescope, index, Stream::pulse)[0];
    jitter = s.pulse_width_ns / fwhm_per_sigma * (u - 0.5);
  }
  p.time_ns = s.pulse_delay_ns - zi * s.radius_m / speed_of_light_m_per_ns + jitter;
  p.wavelength_nm =
      sample(s.spectrum, philox_draw(seed, event, telescope, index, Stream::wavelength));
  p.weight = 1.0;
  return p;
}

inline void emit_background(const BackgroundSource &s, std::uint64_t seed, std::uint32_t event,
                            std::uint32_t telescope, std::uint64_t first, std::size_t rows,
                            double weight, gstrickd::PhotonChunk &chunk) {
  detail::require(s.radius_m > 0.0 && std::isfinite(s.radius_m), "background: radius must be > 0");
  detail::require(s.viewcone_rad > 0.0 && s.viewcone_rad <= pi,
                  "background: viewcone must lie in (0, pi]");
  detail::require(rows <= gstrickd::max_chunk_rows, "a chunk holds at most 65536 rows");
  chunk.resize(rows, source_columns);
  chunk.telescope_id = telescope;
  chunk.frame = gstrickd::Frame::sky;
  chunk.first_index = first;
  using gstrickd::Column;
  for (std::size_t i = 0; i < rows; ++i) {
    const Photon p = sample_background(s, seed, event, telescope, first + i);
    chunk[Column::x][i] = static_cast<float>(p.position.x);
    chunk[Column::y][i] = static_cast<float>(p.position.y);
    chunk[Column::z][i] = static_cast<float>(p.position.z);
    chunk[Column::cx][i] = static_cast<float>(p.direction.x);
    chunk[Column::cy][i] = static_cast<float>(p.direction.y);
    chunk[Column::cz][i] = static_cast<float>(p.direction.z);
    chunk[Column::t][i] = static_cast<float>(p.time_ns);
    chunk[Column::w][i] = static_cast<float>(weight);
    chunk[Column::wavelength][i] = static_cast<float>(p.wavelength_nm);
  }
}

} // namespace dlamb
