#pragma once

// 8-point Gauss-Legendre quadrature shared by the propagation kernels: the node and weight table
// and the piecewise integral with an optional u^2 substitution (the integrands of the near-horizon
// paths have 1/sqrt singularities at the emission or telescope level).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace obdeect::propagation::detail {

// 8-point Gauss-Legendre nodes and weights on [-1, 1].
inline constexpr std::size_t gauss_points = 8;
inline constexpr std::array<double, gauss_points> gauss_nodes{
    -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
    0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
inline constexpr std::array<double, gauss_points> gauss_weights{
    0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
    0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

// Integral over x in [0, length] by Gauss-Legendre on pieces. The first piece has width
// first_width (at most piece_max); each next piece is twice as wide as the previous one, at most
// piece_max. With substitute the pieces are taken in u = sqrt(x) (x = u^2, dx = 2 u du): the
// nodes then resolve a 1/sqrt(x) singularity at x = 0. f(x, weight) is called at each node in
// order and returns false to stop; gauss_pieces then returns false, otherwise true. No allocation.
template <class F>
bool gauss_pieces(double length, double first_width, double piece_max, bool substitute, F &&f) {
  double x0 = 0.;
  double width = std::min(first_width, piece_max);
  while (x0 < length) {
    const double x1 = std::min(length, x0 + width);
    if (substitute) {
      const double u0 = std::sqrt(x0), u1 = std::sqrt(x1);
      const double mid = 0.5 * (u0 + u1), half = 0.5 * (u1 - u0);
      for (std::size_t k = 0; k < gauss_points; ++k) {
        const double u = mid + half * gauss_nodes[k];
        if (!f(u * u, half * 2. * u * gauss_weights[k])) {
          return false;
        }
      }
    } else {
      const double mid = 0.5 * (x0 + x1), half = 0.5 * (x1 - x0);
      for (std::size_t k = 0; k < gauss_points; ++k) {
        if (!f(mid + half * gauss_nodes[k], half * gauss_weights[k])) {
          return false;
        }
      }
    }
    x0 = x1;
    width = std::min(2. * width, piece_max);
  }
  return true;
}

} // namespace obdeect::propagation::detail
