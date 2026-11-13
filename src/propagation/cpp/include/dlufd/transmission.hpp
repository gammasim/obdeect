#pragma once

// Optical-depth table of the atmosphere (extinction-table format) and the transmission of a
// photon from its emission height to the telescope. See docs/extinction.md.

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace dlufd {

// Wavelengths are integer nm in [0, max_lambda).
inline constexpr int max_lambda = 1000;
inline constexpr double unlimited_airmass = 1e10;
inline constexpr double altitude_tolerance_m = 25.0;
inline constexpr double max_telescope_depth = 1000.0;
inline constexpr double max_total_depth = 100.0;
inline constexpr std::size_t max_table_heights = 100;

// Throws unless the wavelength is a finite, non-negative integer.
inline void check_wavelength(double wavelength_nm) {
  if (!std::isfinite(wavelength_nm) || wavelength_nm < 0.0 ||
      wavelength_nm != std::floor(wavelength_nm)) {
    throw std::invalid_argument("wavelength must be a finite non-negative integer (nm)");
  }
}

// Integer wavelength for table lookups, or -1 at or above max_lambda (no transmission).
[[nodiscard]] inline int table_wavelength(double wavelength_nm) {
  check_wavelength(wavelength_nm);
  if (wavelength_nm >= max_lambda) {
    return -1;
  }
  return static_cast<int>(wavelength_nm);
}

namespace detail {

// Text of a line before a '#'.
[[nodiscard]] inline std::string_view trim_comment(std::string_view line) {
  return line.substr(0, line.find('#'));
}

[[nodiscard]] inline std::vector<std::string_view> split_whitespace(std::string_view text) {
  std::vector<std::string_view> tokens;
  std::size_t i = 0;
  while (i < text.size()) {
    while (i < text.size() && std::isspace(static_cast<unsigned char>(text[i]))) {
      ++i;
    }
    const std::size_t start = i;
    while (i < text.size() && !std::isspace(static_cast<unsigned char>(text[i]))) {
      ++i;
    }
    if (i > start) {
      tokens.push_back(text.substr(start, i - start));
    }
  }
  return tokens;
}

// A finite number filling the whole token.
[[nodiscard]] inline double parse_number(std::string_view token, const char *what) {
  const std::string copy(token);
  char *end = nullptr;
  const double value = std::strtod(copy.c_str(), &end);
  if (copy.empty() || end != copy.c_str() + copy.size() || !std::isfinite(value)) {
    throw std::invalid_argument(std::string(what) + " '" + copy + "' is not a finite number");
  }
  return value;
}

// Linear interpolation of y(x) at `at` on ascending x, clamped to the first and last value.
[[nodiscard]] inline double rpol(const std::vector<double> &x, const double *y, double at) {
  if (at <= x.front()) {
    return y[0];
  }
  if (at >= x.back()) {
    return y[x.size() - 1];
  }
  const auto upper = std::upper_bound(x.begin(), x.end(), at);
  const std::size_t i = static_cast<std::size_t>(upper - x.begin());
  return y[i - 1] + (at - x[i - 1]) / (x[i] - x[i - 1]) * (y[i] - y[i - 1]);
}

} // namespace detail

// Optical depth inside the telescope for every integer wavelength.
struct TelescopeDepths {
  std::array<double, max_lambda> od_tel{};   // vertical depth observation level to telescope
  std::array<double, max_lambda> od_focus{}; // vertical depth over the path inside the telescope
};

class TransmissionTable {
public:
  TransmissionTable(const std::string &text, double site_altitude_m) {
    if (!std::isfinite(site_altitude_m)) {
      throw std::invalid_argument("site altitude must be finite");
    }
    parse(text);
    if (std::fabs(1000.0 * h2_km_ - site_altitude_m) > altitude_tolerance_m) {
      throw std::invalid_argument("transmission table H2 does not match the site altitude");
    }
  }

  [[nodiscard]] static TransmissionTable from_file(const std::string &path,
                                                   double site_altitude_m) {
    std::ifstream file(path);
    if (!file) {
      throw std::invalid_argument("cannot open transmission table " + path);
    }
    std::ostringstream text;
    text << file.rdbuf();
    return TransmissionTable(text.str(), site_altitude_m);
  }

  [[nodiscard]] double observation_level_m() const noexcept { return 1000.0 * h2_km_; }
  [[nodiscard]] int first_wavelength_nm() const noexcept { return first_; }
  [[nodiscard]] int last_wavelength_nm() const noexcept { return first_ + rows_ - 1; }

  // Vertical optical depth from the observation level to height_m; linear in log10(height)
  // between the table heights, zero at and below the observation level, constant above the
  // last table height. Wavelengths are truncated to integer nm and clamped to the table.
  template <class T> [[nodiscard]] T optical_depth(T wavelength_nm, T height_m) const {
    const std::size_t row = row_of(static_cast<double>(wavelength_nm));
    const double h = static_cast<double>(height_m);
    if (!std::isfinite(h)) {
      throw std::invalid_argument("height must be finite");
    }
    const double *tau = &tau_[row * grid_.size()];
    if (h <= 1000.0 * h2_km_) {
      return T(0);
    }
    return static_cast<T>(detail::rpol(grid_, tau, std::log10(h)));
  }

  // Transmission exp(-|tau(zem) - od_tel| * airmass - od_focus). Zero below the observation
  // level, for od_tel above 1000 and for a total depth above 100.
  template <class T>
  [[nodiscard]] T transmission2(T wavelength_nm, T zem_m, T airmass, T od_tel, T od_focus) const {
    if (!std::isfinite(static_cast<double>(zem_m))) {
      throw std::invalid_argument("height must be finite");
    }
    if (!std::isfinite(static_cast<double>(airmass)) || airmass < T(0) ||
        !std::isfinite(static_cast<double>(od_tel)) || od_tel < T(0) ||
        !std::isfinite(static_cast<double>(od_focus)) || od_focus < T(0)) {
      throw std::invalid_argument("airmass and optical depths must be finite and non-negative");
    }
    if (static_cast<double>(zem_m) < 1000.0 * h2_km_ ||
        static_cast<double>(od_tel) > max_telescope_depth) {
      (void)row_of(static_cast<double>(wavelength_nm)); // validates the wavelength
      return T(0);
    }
    const T depth = std::fabs(optical_depth<T>(wavelength_nm, zem_m) - od_tel) * airmass + od_focus;
    if (depth > T(max_total_depth)) {
      return T(0);
    }
    return std::exp(-depth);
  }

  template <class T> [[nodiscard]] T transmission(T wavelength_nm, T zem_m, T airmass) const {
    return transmission2<T>(wavelength_nm, zem_m, airmass, T(0), T(0));
  }

  // Depths inside the telescope at altitude telescope_altitude_m with an optical path
  // pathlen_m: od_tel = tau(h_T), od_focus = tau(h_T + l) - tau(h_T). If tau(400 nm, h_T)
  // exceeds 1000 the table does not belong to the site and both are zero.
  [[nodiscard]] TelescopeDepths telescope_depths(double telescope_altitude_m,
                                                 double pathlen_m) const {
    if (!std::isfinite(telescope_altitude_m) || !std::isfinite(pathlen_m) || pathlen_m < 0.0) {
      throw std::invalid_argument("telescope altitude and path length must be finite, path >= 0");
    }
    TelescopeDepths depths;
    if (optical_depth<double>(400.0, telescope_altitude_m) > max_telescope_depth) {
      return depths;
    }
    for (int wavelength = 0; wavelength < max_lambda; ++wavelength) {
      const double w = static_cast<double>(wavelength);
      const double at_telescope = optical_depth<double>(w, telescope_altitude_m);
      depths.od_tel[static_cast<std::size_t>(wavelength)] = at_telescope;
      depths.od_focus[static_cast<std::size_t>(wavelength)] =
          optical_depth<double>(w, telescope_altitude_m + pathlen_m) - at_telescope;
    }
    return depths;
  }

private:
  [[nodiscard]] std::size_t row_of(double wavelength_nm) const {
    if (!std::isfinite(wavelength_nm)) {
      throw std::invalid_argument("wavelength must be finite");
    }
    const double clamped = std::clamp(wavelength_nm, static_cast<double>(first_),
                                      static_cast<double>(last_wavelength_nm()));
    return static_cast<std::size_t>(static_cast<int>(clamped) - first_);
  }

  static double number(const std::string &token, std::size_t line, const char *what) {
    char *end = nullptr;
    const double value = std::strtod(token.c_str(), &end);
    if (token.empty() || end != token.c_str() + token.size() || !std::isfinite(value)) {
      throw std::invalid_argument("transmission table line " + std::to_string(line) + ": " + what +
                                  " '" + token + "' is not a finite number");
    }
    return value;
  }

  void parse(const std::string &text) {
    std::istringstream stream(text);
    std::string line;
    std::size_t line_number = 0;
    std::vector<double> heights_km;
    bool have_header = false;
    std::vector<double> values;
    while (std::getline(stream, line)) {
      ++line_number;
      const std::size_t h2 = line.find("H2=");
      const std::size_t h1 = line.find("H1=");
      if (!have_header && h2 != std::string::npos && h1 != std::string::npos && h2 < h1) {
        std::string head = line.substr(h2 + 3, h1 - h2 - 3);
        std::replace(head.begin(), head.end(), ',', ' ');
        std::istringstream head_stream(head);
        std::string token;
        if (!(head_stream >> token)) {
          throw std::invalid_argument("transmission table: empty H2");
        }
        h2_km_ = number(token, line_number, "H2");
        const std::string_view after_h1 =
            detail::trim_comment(std::string_view(line).substr(h1 + 3));
        std::istringstream list{std::string(after_h1)};
        while (list >> token) {
          heights_km.push_back(number(token, line_number, "height"));
        }
        have_header = true;
        continue;
      }
      std::istringstream row{std::string(detail::trim_comment(line))};
      std::string token;
      if (!(row >> token)) {
        continue;
      }
      if (!have_header) {
        throw std::invalid_argument(
            "transmission table: data before the '# H2= ... H1= ...' header");
      }
      const double wavelength = number(token, line_number, "wavelength");
      if (wavelength != std::floor(wavelength) || wavelength < 0.0 || wavelength >= max_lambda) {
        throw std::invalid_argument("transmission table line " + std::to_string(line_number) +
                                    ": wavelength must be an integer in [0, 1000) nm");
      }
      const int w = static_cast<int>(wavelength);
      if (rows_ == 0) {
        first_ = w;
      } else if (w != first_ + rows_) {
        throw std::invalid_argument("transmission table line " + std::to_string(line_number) +
                                    ": wavelengths must rise in steps of 1 nm");
      }
      const std::size_t before = values.size();
      while (row >> token) {
        const double depth = number(token, line_number, "optical depth");
        if (depth < 0.0) {
          throw std::invalid_argument("transmission table line " + std::to_string(line_number) +
                                      ": optical depth must not be negative");
        }
        values.push_back(depth);
      }
      if (values.size() - before != heights_km.size()) {
        throw std::invalid_argument("transmission table line " + std::to_string(line_number) +
                                    ": number of values differs from the number of heights");
      }
      ++rows_;
    }
    if (heights_km.size() + 1 > max_table_heights) { // H2 counts as a height
      throw std::invalid_argument(
          "transmission table: at most 100 heights (H2 and H1) are supported");
    }
    if (!have_header || heights_km.empty() || rows_ == 0) {
      throw std::invalid_argument("transmission table needs a header with H2 and H1 and data rows");
    }
    if (!(h2_km_ > 0.0)) {
      throw std::invalid_argument("transmission table: H2 must be positive");
    }
    double previous = h2_km_;
    for (double height : heights_km) {
      if (!(height > previous)) {
        throw std::invalid_argument("transmission table: heights must rise above H2");
      }
      previous = height;
    }
    grid_.push_back(std::log10(1000.0 * h2_km_));
    for (double height : heights_km) {
      grid_.push_back(std::log10(1000.0 * height));
    }
    const std::size_t columns = heights_km.size();
    tau_.assign(static_cast<std::size_t>(rows_) * grid_.size(), 0.0);
    for (std::size_t r = 0; r < static_cast<std::size_t>(rows_); ++r) {
      for (std::size_t c = 0; c < columns; ++c) {
        tau_[r * grid_.size() + c + 1] = values[r * columns + c];
      }
    }
  }

  double h2_km_{};
  int first_{};
  int rows_{};
  std::vector<double> grid_; // log10(height / m), starting at the observation level
  std::vector<double> tau_;  // rows x grid_.size(); first column 0
};

} // namespace dlufd
