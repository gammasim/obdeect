#pragma once

// The single cm <-> m conversion used for telescope geometry (TELPOS and CORSIKA positions are in
// cm, the stream in m).

namespace obdeect::stream {

[[nodiscard]] constexpr double cm_to_m(double value_cm) { return value_cm / 100.0; }
[[nodiscard]] constexpr double m_to_cm(double value_m) { return 100.0 * value_m; }

} // namespace obdeect::stream
