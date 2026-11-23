#pragma once

// Optical depth along a straight chord through the atmosphere: the integral of the extinction
// coefficient k(h(s)) = |d tau / dh| over the path length s. Used where the plane formula
// |tau(h_1) - tau(h_2)| / |cos| is ill-conditioned (near-horizontal paths, cos -> 0). See
// docs/propagation/extinction.md, "Chord integral".

#include "obdeect/propagation/transmission.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace obdeect::propagation {

// Paths with |cos| of the local vertical below this use the chord integral instead of the plane
// formula (extinction of sphere rows); scattered photons with |h_s - h_c| / distance below this
// use it in sphere runs.
inline constexpr double near_horizontal_cos = 0.1;

// Largest altitude change [m] of one piece of a chord in the Gauss-Legendre quadrature. The
// extinction coefficient depends on the altitude only, so the pieces are sized by the altitude
// change they cover (the coefficient then changes by a few per cent per piece): a horizontal chord
// is one piece, however long.
inline constexpr double chord_piece_m = 250.0;

// Breakpoints of a chord (path lengths in (0, length), sorted): the coefficient has kinks where the
// chord crosses a table height, and the altitude of a spherical chord is monotone between its
// lowest point and its ends, so the quadrature is split at both. At most two crossings per table
// height, and the lowest point.
inline constexpr std::size_t max_chord_breaks = 2 * max_table_heights + 1;

// Integral of f(s) over s in [0, length_m] (m): 8-point Gauss-Legendre on equal pieces, their
// number max(1, ceil(altitude_change_m / chord_piece_m)) for a chord whose altitude changes
// monotonically by altitude_change_m (absolute value). f is called at the nodes only; no
// allocation. Throws if the length or the altitude change is not finite and non-negative.
template <class F>
[[nodiscard]] double chord_integral(double length_m, double altitude_change_m, F &&f) {
  if (!std::isfinite(length_m) || length_m < 0.0 || !std::isfinite(altitude_change_m) ||
      altitude_change_m < 0.0) {
    throw std::invalid_argument("chord length and altitude change must be finite, not negative");
  }
  if (length_m == 0.0) {
    return 0.0;
  }
  const auto pieces = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::ceil(altitude_change_m / chord_piece_m)));
  const double piece = length_m / static_cast<double>(pieces);
  const double half = 0.5 * piece;
  double sum = 0.0;
  for (std::size_t p = 0; p < pieces; ++p) {
    const double mid = (static_cast<double>(p) + 0.5) * piece;
    for (std::size_t j = 0; j < detail::gauss_points; ++j) {
      sum += detail::gauss_weights[j] * f(mid + half * detail::gauss_nodes[j]);
    }
  }
  return sum * half;
}

// chord_integral over [0, length_m] split at the sorted breaks[0, count) (the breaks in (0,
// length_m); the callers filter them). altitude(s) is the altitude [m] along the chord: within each
// part it must be monotone (the breaks include the extrema); each part is integrated by
// chord_integral with its own altitude change.
template <class Altitude, class F>
[[nodiscard]] double chord_integral_split(double length_m,
                                          const std::array<double, max_chord_breaks> &breaks,
                                          std::size_t count, Altitude &&altitude, F &&f) {
  if (!std::isfinite(length_m) || length_m < 0.0) {
    throw std::invalid_argument("chord length must be finite and not negative");
  }
  double sum = 0.0;
  double start = 0.0;
  for (std::size_t i = 0; i <= count; ++i) {
    const double end = i < count ? breaks[i] : length_m;
    if (end > start) {
      const double change = std::abs(altitude(end) - altitude(start));
      sum += chord_integral(end - start, change, [&](double t) { return f(start + t); });
    }
    start = std::max(start, end);
  }
  return sum;
}

// Flat geometry: tau = integral over the chord of k(h(s)) ds, the height h(s) linear from h_start
// at s = 0 to h_end at s = length_m [m a.s.l.]. For a constant-height chord this is k(h) length.
[[nodiscard]] inline double flat_chord_depth(const TransmissionTable &table, double wavelength_nm,
                                             double h_start_m, double h_end_m, double length_m) {
  if (!std::isfinite(h_start_m) || !std::isfinite(h_end_m)) {
    throw std::invalid_argument("chord heights must be finite");
  }
  if (!(length_m > 0.0)) {
    return 0.0;
  }
  const double slope = (h_end_m - h_start_m) / length_m;
  const auto altitude = [&](double s) { return h_start_m + slope * s; };
  std::array<double, max_chord_breaks> breaks{};
  std::size_t count = 0;
  if (slope != 0.0) {
    for (const double node : table.node_heights_m()) {
      const double s = (node - h_start_m) / slope;
      if (s > 0.0 && s < length_m) {
        breaks[count++] = s;
      }
    }
    std::sort(breaks.begin(), breaks.begin() + static_cast<std::ptrdiff_t>(count));
  }
  return chord_integral_split(length_m, breaks, count, altitude, [&](double s) {
    return table.extinction_coefficient<double>(wavelength_nm, altitude(s));
  });
}

// Spherical geometry (Earth radius earth_radius_m): the straight line from the start point v
// (Earth-centred coordinates [m]) with unit direction d, of length length_m. Its radius is
// r(s)^2 = b^2 + (s + c)^2 with c = v . d and b = |v x d| (the impact parameter, computed from the
// cross product: no cancellation for near-horizontal lines), and its altitude r - R. The chord can
// dip below its ends (the minimum at s = -c); the integral follows the dip. The quadrature is split
// at the lowest point and where the chord crosses a table height. Throws if the length is invalid
// or v is not finite.
[[nodiscard]] inline double spherical_chord_depth(const TransmissionTable &table,
                                                  double wavelength_nm,
                                                  const std::array<double, 3> &v,
                                                  const std::array<double, 3> &d, double length_m,
                                                  double earth_radius_m) {
  for (const double component : v) {
    if (!std::isfinite(component)) {
      throw std::invalid_argument("chord start point must be finite");
    }
  }
  if (!(earth_radius_m > 0.0) || !std::isfinite(earth_radius_m)) {
    throw std::invalid_argument("earth radius must be finite and positive");
  }
  if (!(length_m > 0.0)) {
    return 0.0;
  }
  const double c = v[0] * d[0] + v[1] * d[1] + v[2] * d[2];
  const double cx = v[1] * d[2] - v[2] * d[1];
  const double cy = v[2] * d[0] - v[0] * d[2];
  const double cz = v[0] * d[1] - v[1] * d[0];
  const double b2 = cx * cx + cy * cy + cz * cz;
  const auto altitude = [&](double s) {
    return std::sqrt(b2 + (s + c) * (s + c)) - earth_radius_m;
  };
  std::array<double, max_chord_breaks> breaks{};
  std::size_t count = 0;
  if (-c > 0.0 && -c < length_m) {
    breaks[count++] = -c; // the lowest point
  }
  for (const double node : table.node_heights_m()) {
    const double radius = earth_radius_m + node;
    const double disc = radius * radius - b2;
    if (!(disc >= 0.0)) {
      continue; // the chord does not reach this height
    }
    const double root = std::sqrt(disc);
    for (const double s : {-c - root, -c + root}) {
      if (s > 0.0 && s < length_m) {
        breaks[count++] = s;
      }
    }
  }
  std::sort(breaks.begin(), breaks.begin() + static_cast<std::ptrdiff_t>(count));
  return chord_integral_split(length_m, breaks, count, altitude, [&](double s) {
    return table.extinction_coefficient<double>(wavelength_nm, altitude(s));
  });
}

} // namespace obdeect::propagation
