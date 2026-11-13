#pragma once

// A telescope as listed in a telescope file or in TELPOS: position and radius in cm, and the
// single cm <-> m conversion used for telescope geometry.

namespace obdeect::stream {

struct TelposCm {
  double x{}, y{}, z{}, r{};
};

[[nodiscard]] constexpr double cm_to_m(double value_cm) { return value_cm / 100.0; }
[[nodiscard]] constexpr double m_to_cm(double value_m) { return 100.0 * value_m; }

} // namespace obdeect::stream
