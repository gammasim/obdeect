#pragma once

// Night-sky background source (port of sim_telarray LightEmission BackgroundLightSource,
// used by night-sky). Photons start on the sphere of one telescope, on the hemisphere that
// faces the sky, and travel in a cone of half angle `viewcone_rad` about `pointing`.

#include "dlamb/sky_sources.hpp"
#include "dlamb/source.hpp"
#include "dlamb/streams.hpp"
#include "dlamb/vec3.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/random.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace dlamb {

struct BackgroundSource {
  Vec3 centre{};                           // telescope centre (sky frame, metres)
  double radius_m{};                       // telescope fiducial sphere radius
  Vec3 pointing{0.0, 0.0, -1.0};           // viewing axis (photon directions point away from it)
  double viewcone_rad{1.5707963267948966}; // half opening angle of the cone
  Spectrum spectrum{Fixed{400.0}};
  double pulse_width_ns{0.0}; // width of the flat time distribution (as original)
  double pulse_delay_ns{0.0};
};

// Mean of bunches for one event and telescope: floor plus one with the fractional
// probability (LightEmission bunch rounding). Keys as photon counts of dlamb/source.hpp.
[[nodiscard]] inline std::uint64_t background_count(double mean_bunches, std::uint64_t seed,
                                                    std::uint32_t event, std::uint32_t telescope) {
  const double whole = std::floor(mean_bunches);
  const auto u = gstrickd::uniform4(
      {seed, stage_sources, event, telescope, 0, static_cast<std::uint32_t>(Stream::count)});
  return static_cast<std::uint64_t>(whole) + (u[0] < mean_bunches - whole ? 1U : 0U);
}

// Photon number `index` of one background source. Directions are drawn uniformly in
// cos(theta) inside the cone; the emission point follows the LightEmission construction
// (projected radius 1 - sqrt(U) on the unit hemisphere, rotated with Rotate_y, Rotate_z),
// and its time is shifted by the travel time from the sphere to the observation level.
[[nodiscard]] inline Photon sample_background(const BackgroundSource &s, std::uint64_t seed,
                                              std::uint32_t event, std::uint32_t telescope,
                                              std::uint64_t index) {
  const auto key = [&](Stream stream) {
    return gstrickd::DrawKey{seed,      stage_sources, event,
                             telescope, index,         static_cast<std::uint32_t>(stream)};
  };
  Photon photon;
  const Uniforms du = gstrickd::uniform4(key(Stream::direction));
  const double cos_theta = 1.0 - du[0] * (1.0 - std::cos(s.viewcone_rad));
  const Vec3 d = detail::around_axis(normalised(s.pointing), cos_theta, two_pi * du[1]);
  photon.direction = d;

  const Uniforms tu = gstrickd::uniform4(key(Stream::target));
  const double ri = 1.0 - std::sqrt(tu[0]);
  const double phi = two_pi * tu[1];
  const double zi = std::sqrt(std::max(0.0, 1.0 - ri * ri));
  // The hemisphere is centred on -d (the direction towards the source), so every point
  // u on it has dot(u, d) <= 0 and the photon enters the sphere. Both angles belong to -d.
  const Vec3 towards_source = -d;
  const double za = std::acos(std::clamp(towards_source.z, -1.0, 1.0));
  const double az = std::atan2(towards_source.y, towards_source.x);
  const Vec3 unit = rotate_z(rotate_y({ri * std::cos(phi), ri * std::sin(phi), zi}, za), az);
  photon.position = s.centre + unit * s.radius_m;

  const double travel_ns = zi * s.radius_m / speed_of_light_m_per_ns;
  const Uniforms pu = gstrickd::uniform4(key(Stream::pulse));
  const double width =
      s.pulse_width_ns > 0.0 ? s.pulse_width_ns / 2.35482004503095 * (pu[0] - 0.5) : 0.0;
  photon.time_ns = s.pulse_delay_ns - travel_ns + width;

  const Uniforms wu = gstrickd::uniform4(key(Stream::wavelength));
  photon.wavelength_nm = std::visit([&](const auto &sp) { return sp.sample(wu); }, s.spectrum);
  return photon;
}

// Rows first_index .. first_index + rows - 1 of one background source, with columns of
// dlamb::source_columns; every row has weight `weight` (the bunch size).
inline void emit_background(const BackgroundSource &s, std::uint64_t seed, std::uint32_t event,
                            std::uint32_t telescope, std::uint64_t first_index, std::size_t rows,
                            double weight, gstrickd::PhotonChunk &out) {
  if (!(s.radius_m > 0.0) || !(s.viewcone_rad > 0.0) || s.viewcone_rad > pi)
    throw std::invalid_argument("background source needs radius > 0 and a cone in (0, 180] deg");
  if (rows > gstrickd::max_chunk_rows)
    throw std::invalid_argument("a chunk holds at most 2^16 photons");
  out.telescope_id = telescope;
  out.frame = gstrickd::Frame::sky;
  out.first_index = first_index;
  out.resize(rows, source_columns);
  auto &x = out[gstrickd::Column::x], &y = out[gstrickd::Column::y], &z = out[gstrickd::Column::z];
  auto &cx = out[gstrickd::Column::cx], &cy = out[gstrickd::Column::cy],
       &cz = out[gstrickd::Column::cz];
  auto &t = out[gstrickd::Column::t], &w = out[gstrickd::Column::w],
       &wl = out[gstrickd::Column::wavelength];
  for (std::size_t i = 0; i < rows; ++i) {
    const Photon p = sample_background(s, seed, event, telescope, first_index + i);
    x[i] = static_cast<float>(p.position.x);
    y[i] = static_cast<float>(p.position.y);
    z[i] = static_cast<float>(p.position.z);
    cx[i] = static_cast<float>(p.direction.x);
    cy[i] = static_cast<float>(p.direction.y);
    cz[i] = static_cast<float>(p.direction.z);
    t[i] = static_cast<float>(p.time_ns);
    w[i] = static_cast<float>(weight);
    wl[i] = static_cast<float>(p.wavelength_nm);
  }
}

} // namespace dlamb
