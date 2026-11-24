#pragma once

#include "obdeect/ray_tracing/model.hpp"

#include <algorithm>
#include <cstddef>

namespace obdeect::ray_tracing {

namespace detail {

// Upper bound of Response::at(wavelength, any angle): the angle interpolation
// is a convex combination of node values, so the largest node of the two
// interpolated wavelength rows bounds it. Same wavelength arithmetic as
// Response::at.
inline double response_bound(const Response &r, double wavelength_nm) {
  const double position = std::clamp((wavelength_nm - r.wavelength_min_nm) / r.wavelength_step_nm,
                                     0.0, static_cast<double>(r.wavelength_count - 1));
  const auto w = std::min(static_cast<std::size_t>(position), r.wavelength_count - 1);
  const double fw = position - static_cast<double>(w);
  const std::size_t w1 = std::min(w + 1, r.wavelength_count - 1);
  const std::size_t na = r.incidence_deg.size();
  double bound = r.values[w * na] * (1 - fw) + r.values[w1 * na] * fw;
  for (std::size_t j = 1; j < na; ++j)
    bound = std::max(bound, r.values[w * na + j] * (1 - fw) + r.values[w1 * na + j] * fw);
  return bound;
}

// Largest node of a degradation map (the factor of its bound); 1 without a map.
inline double degradation_maximum(const Mirror &mirror) {
  return mirror.degradation ? *std::max_element(mirror.degradation->values.begin(),
                                                mirror.degradation->values.end())
                            : 1.0;
}

// Bound of the mirror's response at a wavelength, given degradation_maximum(mirror).
// Bilinear interpolation never exceeds the largest node; outside the map the value is
// either the edge value or 0.
inline double response_bound(const Mirror &mirror, double wavelength_nm, double degradation_max) {
  double bound = 1.0;
  if (mirror.reflectivity)
    bound *= response_bound(*mirror.reflectivity, wavelength_nm);
  if (mirror.degradation)
    bound *= degradation_max;
  return bound;
}

} // namespace detail

// Throws ModelError for a telescope transmission whose maximum is not known.
inline void validate_bound_support(const OpticalModel &model) {
  if (const auto &t = model.telescope_transmission)
    if (!(t->amplitude >= 0.0 && t->power > 0.0 && t->angular_scale_rad > 0.0 &&
          t->outer_power >= 0.0))
      throw ModelError("response_bound: unsupported telescope transmission");
}

// An upper bound of the weight factor that trace_photon with apply_responses
// can apply to a photon of this wavelength: reflectivities, degradation maps
// and the telescope transmission. Throws ModelError for a transmission whose
// maximum is not known.
//
// The bound is the scale of the bunch roulette: a photon is kept with probability
// weight / bound (at most one). A photon whose weight exceeds the bound (a direct
// detection without reflection, or a detection that reflects on fewer mirrors than
// the sequence) is kept with probability one and carries the excess in its weight,
// so the estimator stays unbiased (trace_bunches). The bound is tight for the
// sequential photons, which therefore keep their output.
// The same bound for many wavelengths of one model: the degradation maxima are computed once.
// The model must outlive the object.
class ResponseBound {
public:
  explicit ResponseBound(const OpticalModel &model)
      : model_(model), primary_degradation_(detail::degradation_maximum(model.primary)),
        secondary_degradation_(model.secondary ? detail::degradation_maximum(*model.secondary)
                                               : 1.0) {
    validate_bound_support(model);
  }

  [[nodiscard]] double at(double wavelength_nm) const {
    double bound = detail::response_bound(model_.primary, wavelength_nm, primary_degradation_);
    if (model_.secondary)
      bound *= detail::response_bound(*model_.secondary, wavelength_nm, secondary_degradation_);
    if (const auto &t = model_.telescope_transmission)
      bound *= t->on_axis;
    if (model_.camera_window)
      bound *= detail::response_bound(model_.camera_window->internal_transmittance, wavelength_nm);
    return bound * (1.0 + 1e-12);
  }

private:
  const OpticalModel &model_;
  double primary_degradation_, secondary_degradation_;
};

[[nodiscard]] inline double response_bound(const OpticalModel &model, double wavelength_nm) {
  return ResponseBound(model).at(wavelength_nm);
}

} // namespace obdeect::ray_tracing
