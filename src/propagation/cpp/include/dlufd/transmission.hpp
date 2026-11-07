#pragma once

// Atmospheric transmission tables and telescope optical depths, as in sim_telarray
// (sim_absorb.c: read_trans, atmospheric_transmission, atmospheric_transmission2;
// rpolator.c: rpol; sim_telarray.c: od2tel and od2focus, MAX_LAMBDA 1000).
//
// Units: altitudes in metres a.s.l. at the API; sim_telarray works in cm, and
// the conversion is a multiplication by 100 before the table lookup.

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace dlufd {

// Limits and sentinels of sim_telarray.
inline constexpr int max_lambda_nm = 1000;                 // MAX_LAMBDA (mc_aux.h)
inline constexpr std::size_t max_table_heights = 100;      // MAX_TRANS_HEIGHT (sim_absorb.c)
inline constexpr std::size_t max_table_wavelengths = 1200; // MAX_WL (sim_absorb.c)
inline constexpr double optical_depth_limit = 100.0;       // transmission cut (sim_absorb.c)
inline constexpr double telescope_optical_depth_limit = 1000.0;
inline constexpr double unset_optical_depth = 99999.0; // sim_telarray.c: od when t <= 0
inline constexpr double altitude_tolerance_m = 25.0;   // read_trans altitude check
inline constexpr double unlimited_airmass = 1e10;      // sim_telarray.c: cz == 0

// Optical depths per integer wavelength (0 .. max_lambda_nm - 1) for one telescope:
// od_tel between the telescope and the observation level, od_focus inside the telescope.
struct TelescopeDepths {
  std::array<double, max_lambda_nm> od_tel{};
  std::array<double, max_lambda_nm> od_focus{};
};

namespace detail {

inline std::string_view trim_comment(std::string_view line) {
  const auto hash = line.find('#');
  return hash == std::string_view::npos ? line : line.substr(0, hash);
}

inline std::vector<std::string_view> split_whitespace(std::string_view text) {
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

// C atoi(): optional sign and leading digits, 0 if there are none.
inline long c_atoi(std::string_view token) {
  std::size_t i = 0;
  bool negative = false;
  if (i < token.size() && (token[i] == '+' || token[i] == '-')) {
    negative = token[i] == '-';
    ++i;
  }
  long value = 0;
  for (; i < token.size() && std::isdigit(static_cast<unsigned char>(token[i])); ++i) {
    value = value * 10 + (token[i] - '0');
  }
  return negative ? -value : value;
}

inline double parse_number(std::string_view token, const char *what) {
  const std::string text(token);
  char *end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  if (text.empty() || end != text.c_str() + text.size()) {
    throw std::invalid_argument(std::string("invalid number in ") + what + ": '" + text + "'");
  }
  return value;
}

// rpol() of rpolator.c: linear interpolation in ascending x, clamped to the edge
// values outside [x[0], x[n-1]]. The result for an interior point does not depend on
// which of two intervals sharing a node is used, so a binary search gives the C result.
template <std::floating_point T>
[[nodiscard]] T rpol(const std::vector<double> &x, const double *y, T xp) {
  const std::size_t n = x.size();
  if (static_cast<double>(xp) <= x.front()) {
    return T{y[0]};
  }
  if (static_cast<double>(xp) >= x.back()) {
    return T{y[n - 1]};
  }
  const auto hi = static_cast<std::size_t>(
      std::upper_bound(x.begin(), x.end(), static_cast<double>(xp)) - x.begin());
  const std::size_t lo = hi - 1;
  const T rpl = (xp - T{x[lo]}) / T{x[hi] - x[lo]};
  return T{y[lo]} * (T{1} - rpl) + T{y[hi]} * rpl;
}

} // namespace detail

// A sim_telarray atmospheric transmission table: optical depth as a function of
// wavelength (integer nm rows) and emission altitude (log10 of height in cm).
class TransmissionTable {
public:
  // Parse the text of a sim_telarray transmission file (read_trans, sim_absorb.c:84).
  // site_altitude_m is the observation level; it must match the table's H2 within 25 m.
  // Throws std::invalid_argument on invalid input.
  TransmissionTable(std::string_view text, double site_altitude_m) { parse(text, site_altitude_m); }

  // Read a transmission file from disk. Throws std::runtime_error if it cannot be read.
  static TransmissionTable from_file(const std::string &path, double site_altitude_m) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
      throw std::runtime_error("cannot read transmission file: " + path);
    }
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return TransmissionTable(buffer.str(), site_altitude_m);
  }

  [[nodiscard]] double observation_level_m() const noexcept { return h2_km_ * 1000.0; }
  [[nodiscard]] int first_wavelength_nm() const noexcept { return wl_first_; }
  [[nodiscard]] int last_wavelength_nm() const noexcept {
    return wl_first_ + static_cast<int>(rows_.size()) - 1;
  }

  // Optical depth of a vertical path from the observation level to height_m, at
  // the wavelength (truncated to integer nm and clamped to the table range).
  template <std::floating_point T>
  [[nodiscard]] T optical_depth(double lambda_nm, T height_m) const {
    return detail::rpol(log_h_cm_, row(lambda_nm).data(), std::log10(height_m * T{100}));
  }

  // sim_telarray atmospheric_transmission(): exp(-od * airmass) for an emission at
  // zem_m a.s.l., or 0 below the observation level or above the optical-depth cut.
  template <std::floating_point T>
  [[nodiscard]] T transmission(double lambda_nm, T zem_m, T airmass) const {
    const T zem_cm = zem_m * T{100};
    if (zem_cm < T{h2_cm()} || zem_cm <= T{0}) {
      return T{0};
    }
    const T od = detail::rpol(log_h_cm_, row(lambda_nm).data(), std::log10(zem_cm)) * airmass;
    if (od > T{optical_depth_limit}) {
      return T{0};
    }
    return std::exp(-od);
  }

  // sim_telarray atmospheric_transmission2(): the default extinction of photon bunches.
  // od_tel and od_focus are the telescope optical depths at this wavelength
  // (TelescopeDepths); od_focus is not multiplied by the airmass.
  template <std::floating_point T>
  [[nodiscard]] T transmission2(double lambda_nm, T zem_m, T airmass, T od_tel, T od_focus) const {
    const T zem_cm = zem_m * T{100};
    if (zem_cm < T{h2_cm()} || zem_cm <= T{0} || od_tel > T{telescope_optical_depth_limit}) {
      return T{0};
    }
    const T od =
        std::abs(detail::rpol(log_h_cm_, row(lambda_nm).data(), std::log10(zem_cm)) - od_tel) *
            airmass +
        od_focus;
    if (od > T{optical_depth_limit}) {
      return T{0};
    }
    return std::exp(-od);
  }

  // Per-telescope optical depths for all integer wavelengths 0..999, as set up in
  // sim_telarray.c (lines 1854-1890). telescope_altitude_m is the telescope's altitude
  // a.s.l. and pathlen_m the optical path inside the telescope (sim_config.c, pathlen_tel).
  [[nodiscard]] TelescopeDepths telescope_depths(double telescope_altitude_m,
                                                 double pathlen_m) const {
    TelescopeDepths depths;
    for (int iwl = 0; iwl < max_lambda_nm; ++iwl) {
      const double t1 = transmission<double>(iwl, telescope_altitude_m, 1.0);
      const double t2 = transmission<double>(iwl, telescope_altitude_m + pathlen_m, 1.0);
      depths.od_tel[iwl] = t1 <= 0.0 ? unset_optical_depth : -std::log(t1);
      depths.od_focus[iwl] = (t1 <= 0.0 || t2 <= 0.0) ? unset_optical_depth : -std::log(t2 / t1);
    }
    // sim_telarray.c:1866: a 400 nm depth above the limit means the table does not fit
    // the site; sim_telarray then ignores the telescope depths entirely.
    if (depths.od_tel[400] > telescope_optical_depth_limit) {
      depths.od_tel.fill(0.0);
      depths.od_focus.fill(0.0);
    }
    return depths;
  }

private:
  double h2_km_{};
  int wl_first_{};
  std::vector<double> log_h_cm_{};
  std::vector<std::vector<double>> rows_{}; // rows_[iwl - wl_first_][height index]

  [[nodiscard]] double h2_cm() const noexcept { return h2_km_ * 1e5; }

  // Table row for an integer wavelength, clamped to the table range (sim_absorb.c).
  [[nodiscard]] const std::vector<double> &row(double lambda_nm) const {
    const auto iwl = static_cast<long>(lambda_nm);
    const long first = wl_first_;
    const long last = last_wavelength_nm();
    return rows_[static_cast<std::size_t>(std::clamp(iwl, first, last) - first)];
  }

  void parse(std::string_view text, double site_altitude_m) {
    // Header: first line starting with "# H2=".
    std::string_view header;
    std::string_view body;
    bool found = false;
    std::size_t pos = 0;
    while (pos < text.size()) {
      std::size_t end = text.find('\n', pos);
      if (end == std::string_view::npos) {
        end = text.size();
      }
      std::string_view line = text.substr(pos, end - pos);
      if (!line.empty() && line.back() == '\r') {
        line.remove_suffix(1);
      }
      pos = std::min(end + 1, text.size());
      if (line.starts_with("# H2=")) {
        header = line;
        body = text.substr(pos);
        found = true;
        break;
      }
    }
    if (!found) {
      throw std::invalid_argument("transmission table: no '# H2=' header line");
    }

    // H2 is the first number after "# H2=" (sscanf "%lf"; a comma ends it).
    std::string_view h2_field = header.substr(5);
    h2_field = h2_field.substr(0, h2_field.find(','));
    const auto h2_tokens = detail::split_whitespace(h2_field);
    if (h2_tokens.empty()) {
      throw std::invalid_argument("transmission table: missing H2 value");
    }
    h2_km_ = detail::parse_number(h2_tokens[0], "H2");
    if (std::abs(h2_km_ * 1000.0 - site_altitude_m) > altitude_tolerance_m) {
      throw std::invalid_argument("transmission table H2 does not match the site altitude");
    }

    // Heights: H2 followed by the H1 list (sim_absorb.c:198-227).
    const auto h1_at = header.find(" H1=");
    if (h1_at == std::string_view::npos) {
      throw std::invalid_argument("transmission table: missing H1 list");
    }
    std::vector<double> heights_km{h2_km_};
    for (const auto token :
         detail::split_whitespace(detail::trim_comment(header.substr(h1_at + 4)))) {
      heights_km.push_back(detail::parse_number(token, "H1"));
    }
    if (heights_km.size() > max_table_heights) {
      throw std::invalid_argument("transmission table has more heights than MAX_TRANS_HEIGHT");
    }
    if (heights_km.size() < 2) {
      throw std::invalid_argument("transmission table needs at least one H1 height");
    }
    for (std::size_t i = 1; i < heights_km.size(); ++i) {
      if (heights_km[i] <= heights_km[i - 1]) {
        throw std::invalid_argument("transmission table heights are not increasing");
      }
    }
    log_h_cm_.clear();
    for (const double h : heights_km) {
      log_h_cm_.push_back(std::log10(h * 1e5));
    }

    // Rows: wavelength then one optical depth per H1 height. Lines without a wavelength
    // or with too few values are skipped, as in read_trans. Column 0 (H2) is 0.
    std::vector<long> wavelengths;
    std::size_t pos_body = 0;
    while (pos_body < body.size()) {
      std::size_t end = body.find('\n', pos_body);
      if (end == std::string_view::npos) {
        end = body.size();
      }
      const std::string_view line = body.substr(pos_body, end - pos_body);
      pos_body = std::min(end + 1, body.size());
      const auto tokens = detail::split_whitespace(detail::trim_comment(line));
      if (tokens.empty()) {
        continue;
      }
      const long wavelength = detail::c_atoi(tokens[0]);
      if (wavelength == 0 || tokens.size() < heights_km.size()) {
        continue;
      }
      std::vector<double> values(heights_km.size(), 0.0);
      for (std::size_t i = 1; i < heights_km.size(); ++i) {
        values[i] = detail::parse_number(tokens[i], "optical depth");
      }
      wavelengths.push_back(wavelength);
      rows_.push_back(std::move(values));
      if (rows_.size() >= max_table_wavelengths) {
        throw std::invalid_argument("transmission table has MAX_WL or more wavelengths");
      }
    }
    if (rows_.empty()) {
      throw std::invalid_argument("transmission table: no transmission data");
    }
    for (std::size_t i = 0; i < wavelengths.size(); ++i) {
      if (wavelengths[i] != wavelengths[0] + static_cast<long>(i)) {
        throw std::invalid_argument("transmission table: wavelengths not at 1 nm intervals");
      }
    }
    wl_first_ = static_cast<int>(wavelengths[0]);
  }
};

} // namespace dlufd
