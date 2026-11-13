#pragma once

// Cherenkov light of one muon-like track filtered to a telescope sphere, sky frame
// (docs/spec/tool-muon-track.md).

#include "obdeect/sources/cli.hpp"
#include "obdeect/sources/emission_stream.hpp"
#include "obdeect/sources/track.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect::sources::muon_track {

inline constexpr double alpha_fine = 0.0072973526;

// Refractivity n - 1 against altitude, linearly interpolated, constant outside the nodes.
struct IndexProfile {
  std::vector<double> altitude_m;
  std::vector<double> index_minus_one;

  [[nodiscard]] bool empty() const { return altitude_m.empty() && index_minus_one.empty(); }

  [[nodiscard]] double at(double h) const {
    if (h <= altitude_m.front())
      return index_minus_one.front();
    if (h >= altitude_m.back())
      return index_minus_one.back();
    const auto upper = std::upper_bound(altitude_m.begin(), altitude_m.end(), h);
    const std::size_t i = static_cast<std::size_t>(upper - altitude_m.begin()) - 1;
    const double f = (h - altitude_m[i]) / (altitude_m[i + 1] - altitude_m[i]);
    return index_minus_one[i] + f * (index_minus_one[i + 1] - index_minus_one[i]);
  }
};

struct Options {
  std::int64_t run_number{1};
  std::uint64_t seed{1};
  std::uint32_t events{1};
  double observation_level_m{2150.0};
  double telescope_radius_m{15.0};
  double impact_x_m{0.0};
  double impact_y_m{0.0};
  double rscat_m{0.0};
  int nscat{1};
  double h_start_m{3000.0};
  double theta_deg{0.0};
  double phi_deg{0.0};
  double index_minus_one{0.000292};
  IndexProfile index_profile;
  double bunch_size{5.0};
  double wavelength_min_nm{240.0};
  double wavelength_max_nm{1000.0};
  double charge{-1.0};
  double beta{1.0};
  std::optional<double> momentum_gev;
  double mass_gev{0.105658};
  int primary_id{6};
  bool multiple_scattering{false};
};

// Number of events from a command-line integer: 1..4294967295.
inline std::uint32_t event_count(long long n) {
  if (n < 1 || n > 4294967295LL)
    throw std::invalid_argument("events must be between 1 and 4294967295");
  return static_cast<std::uint32_t>(n);
}

inline int primary_id_from_name(const std::string &name) {
  if (name == "gamma")
    return 1;
  if (name == "positron")
    return 2;
  if (name == "electron")
    return 3;
  if (name == "muon" || name == "mu-")
    return 6;
  if (name == "mu+")
    return 5;
  if (name == "proton")
    return 14;
  throw std::invalid_argument("unknown primary name: " + name);
}

// Mass [GeV] and charge of a CORSIKA particle id.
inline void mass_charge(int id, double &mass, double &charge) {
  struct Entry {
    int id;
    double charge, mass;
  };
  static constexpr Entry table[] = {
      {3, -1, 0.000511}, {2, 1, 0.000511},   {6, -1, 0.105658},  {5, 1, 0.105658},
      {8, 1, 0.13957},   {11, 1, 0.493677},  {14, 1, 0.938272},  {13, 0, 0.939565},
      {25, 0, 0.939565}, {9, -1, 0.13957},   {12, -1, 0.493677}, {10, 0, 0.497648},
      {16, 0, 0.497648}, {15, -1, 0.938272}, {66, 0, 1e-11},     {67, 0, 1e-11},
      {68, 0, 1e-11},    {69, 0, 1e-11},     {133, 0, 1e-11},    {134, 0, 1e-11},
      {130, 1, 1.777},   {131, -1, 1.777},   {7, 0, 0.134},      {17, 0, 0.54751},
      {50, 0, 0.78265},  {51, 0, 0.7755},    {52, 1, 0.7755},    {53, -1, 0.7755}};
  for (const auto &e : table)
    if (e.id == id) {
      charge = e.charge;
      mass = e.mass;
      return;
    }
  if (id <= 1) {
    charge = 0.0;
    mass = 0.0;
  } else if ((id > 1 && id < 100) || (id > 110 && id < 199)) {
    charge = 0.0;
    mass = 1.2;
  } else if (id == 402) {
    charge = 2.0;
    mass = 0.931845 * 4.0;
  } else {
    charge = id % 100;
    mass = (0.931494 - 0.000511) * (id / 100);
  }
}

// Whitespace separated "altitude n-1" pairs.
inline IndexProfile read_index_profile(std::istream &in, const std::string &name) {
  std::vector<double> values;
  std::string word;
  while (in >> word) {
    if (!cli::is_number(word))
      throw std::invalid_argument("index profile " + name + ": expected 'altitude n-1' lines");
    values.push_back(std::strtod(word.c_str(), nullptr));
  }
  if (values.empty())
    throw std::invalid_argument("index profile " + name + " has no data");
  if (values.size() % 2 != 0)
    throw std::invalid_argument("index profile " + name + ": a value is missing its pair");
  IndexProfile profile;
  for (std::size_t i = 0; i < values.size(); i += 2) {
    profile.altitude_m.push_back(values[i]);
    profile.index_minus_one.push_back(values[i + 1]);
  }
  return profile;
}

inline void validate(const Options &o) {
  using detail::require;
  require(o.events >= 1, "events must be between 1 and 4294967295");
  require(within_limits(o.observation_level_m), "observation level out of range");
  require(positive_within_limits(o.telescope_radius_m), "telescope radius out of range");
  require(within_limits(o.impact_x_m) && within_limits(o.impact_y_m),
          "impact position out of range");
  require(positive_within_limits(o.h_start_m), "start altitude out of range");
  require(std::isfinite(o.theta_deg) && o.theta_deg >= 0.0 && o.theta_deg < 90.0,
          "theta must lie in [0, 90) deg");
  require(std::isfinite(o.phi_deg), "phi must be finite");
  require(std::isfinite(o.index_minus_one) && o.index_minus_one > 0.0 && o.index_minus_one < 0.0003,
          "n-1 must lie in (0, 0.0003)");
  require(std::isfinite(o.bunch_size) && o.bunch_size > 0.0, "bunch size must be > 0");
  require(std::isfinite(o.wavelength_min_nm) && std::isfinite(o.wavelength_max_nm) &&
              o.wavelength_min_nm >= 200.0 && o.wavelength_min_nm <= o.wavelength_max_nm &&
              o.wavelength_max_nm <= 2000.0,
          "wavelengths must satisfy 200 <= wlmin <= wlmax <= 2000 nm");
  require(std::isfinite(o.charge) && o.charge != 0.0, "charge must be finite and non-zero");
  require(std::isfinite(o.mass_gev) && o.mass_gev >= 0.0005 && o.mass_gev <= 100.0,
          "mass must lie in [0.0005, 100] GeV");
  require(std::isfinite(o.beta) && o.beta > 0.0 && o.beta <= 1.0, "beta must lie in (0, 1]");
  require(!o.momentum_gev || (std::isfinite(*o.momentum_gev) && *o.momentum_gev > 0.0),
          "momentum must be > 0");
  require(!o.multiple_scattering, "multiple scattering is not supported");
  require(o.rscat_m == 0.0, "--rscat is not supported");
  require(o.nscat == 1, "--nscat must be 1 (each event is used once)");
  const auto &p = o.index_profile;
  if (!p.empty()) {
    require(p.altitude_m.size() == p.index_minus_one.size() && p.altitude_m.size() >= 2,
            "index profile needs at least two nodes with equal numbers of altitudes and values");
    for (std::size_t i = 0; i < p.altitude_m.size(); ++i) {
      require(std::isfinite(p.altitude_m[i]) && std::isfinite(p.index_minus_one[i]) &&
                  p.index_minus_one[i] > 0.0,
              "index profile values must be finite, n-1 positive");
      require(i == 0 || p.altitude_m[i] > p.altitude_m[i - 1],
              "index profile altitudes must increase");
    }
  }
}

namespace detail {

inline Track make_track(const Options &o, const Vec3 &start, const Vec3 &direction, double length,
                        double cos_ct, double time_ns) {
  Track t;
  t.start = start;
  t.direction = direction;
  t.length_m = length;
  t.cos_theta = cos_ct;
  t.beta = o.beta;
  t.time_ns = time_ns;
  if (o.wavelength_min_nm == o.wavelength_max_nm)
    t.spectrum = Fixed{o.wavelength_min_nm};
  else
    t.spectrum = InverseSquare{o.wavelength_min_nm, o.wavelength_max_nm};
  t.target = SphereTarget{{0.0, 0.0, 0.0}, o.telescope_radius_m};
  t.photons = 2.0 * pi * alpha_fine * o.charge * o.charge * (1.0 - cos_ct * cos_ct) *
              (1.0 / o.wavelength_min_nm - 1.0 / o.wavelength_max_nm) * 1e9 * length;
  return t;
}

} // namespace detail

// Track segments of the muon, lowest first.
inline std::vector<Track> tracks(const Options &o) {
  validate(o);
  constexpr double deg = pi / 180.0;
  const double theta = o.theta_deg * deg, phi = o.phi_deg * deg;
  const Vec3 direction{-std::sin(theta) * std::cos(phi), -std::sin(theta) * std::sin(phi),
                       -std::cos(theta)};
  const double dist = o.h_start_m / std::cos(theta);
  const Vec3 impact{o.impact_x_m, o.impact_y_m, 0.0};
  std::vector<Track> out;
  if (o.index_profile.empty()) {
    const double cos_ct = 1.0 / ((1.0 + o.index_minus_one) * o.beta);
    if (cos_ct <= 1.0)
      out.push_back(detail::make_track(o, impact - direction * dist, direction, dist, cos_ct,
                                       -dist / (o.beta * speed_of_light_m_per_ns)));
    return out;
  }
  const auto steps = static_cast<std::size_t>(std::ceil(dist / 10.0));
  const double step = dist / static_cast<double>(steps);
  Vec3 lower = impact;
  double travelled = 0.0;
  for (std::size_t k = 0; k < steps; ++k) {
    const Vec3 upper = lower - direction * step;
    const Vec3 mid = lower - direction * (0.5 * step);
    travelled += step;
    const double cos_ct = 1.0 / ((1.0 + o.index_profile.at(mid.z)) * o.beta);
    if (cos_ct <= 1.0)
      out.push_back(detail::make_track(o, upper, direction, step, cos_ct,
                                       -travelled / (o.beta * speed_of_light_m_per_ns)));
    lower = upper;
  }
  return out;
}

inline std::string usage() {
  return "usage: obdeect sources muon-track [options]\n"
         "  -r R  telescope radius [m]      -x X -y Y  impact position [m]\n"
         "  -h H  start altitude [m]        --theta T  zenith angle [deg]   --phi P  azimuth "
         "[deg]\n"
         "  -n V  refractive index (n or n-1)   --index-profile FILE  altitude [m] and n-1\n"
         "  -a M  observation level [m]     -l/-L  wavelength range [nm]\n"
         "  -q Q  charge    -m M  mass [GeV]    -p P  momentum [GeV]    --beta B\n"
         "  --primary ID|NAME   -N N  events   --run N   --seed N   -o FILE   --help\n";
}

struct Parsed {
  Options options;
  std::string output;
  bool help{false};
};

namespace detail {

inline double number(const std::string &text, const std::string &option) {
  if (!cli::is_number(text))
    throw std::invalid_argument("invalid number for " + option + ": '" + text + "'");
  const double v = std::strtod(text.c_str(), nullptr);
  if (!std::isfinite(v))
    throw std::invalid_argument("invalid number for " + option + ": '" + text + "'");
  return v;
}

inline long long integer(const std::string &text, const std::string &option) {
  try {
    return cli::parse_integer(text, option);
  } catch (const std::invalid_argument &) {
    throw std::invalid_argument("invalid integer for " + option + ": '" + text + "'");
  }
}

} // namespace detail

// Parses the arguments (without argv[0]). Throws std::invalid_argument.
inline Parsed parse_command_line(const std::vector<std::string> &args) {
  Parsed parsed;
  Options &o = parsed.options;
  std::optional<double> beta_given;
  for (std::size_t i = 0; i < args.size(); ++i) {
    const std::string &arg = args[i];
    if (arg == "--help") {
      parsed.help = true;
      return parsed;
    }
    if (arg == "--enable-scatter" || arg == "--scatter") {
      o.multiple_scattering = true;
      continue;
    }
    if (arg == "--disable-scatter" || arg == "--no-scatter") {
      o.multiple_scattering = false;
      continue;
    }
    if (arg.rfind("-I", 0) == 0)
      throw std::invalid_argument("-I (atmosphere table search path) is not supported");
    const auto is_known = [&] {
      static const char *const names[] = {"-r",
                                          "--radius",
                                          "--telescope-radius",
                                          "-x",
                                          "-y",
                                          "-h",
                                          "--h-start",
                                          "--theta",
                                          "--zenith-angle",
                                          "--phi",
                                          "--orientation",
                                          "--azimuth",
                                          "-n",
                                          "--n1",
                                          "--refraction",
                                          "-A",
                                          "--atmosphere",
                                          "--index-profile",
                                          "-a",
                                          "--altitude",
                                          "--obslev",
                                          "--observation-level",
                                          "-b",
                                          "--bs",
                                          "--bunch-size",
                                          "-l",
                                          "--wlmin",
                                          "-L",
                                          "--wlmax",
                                          "-q",
                                          "--charge",
                                          "-m",
                                          "--mass",
                                          "-p",
                                          "--momentum",
                                          "--beta",
                                          "--primary",
                                          "-N",
                                          "--events",
                                          "--run",
                                          "--seed",
                                          "-o",
                                          "--output",
                                          "--rscat",
                                          "--nscat"};
      return std::find(std::begin(names), std::end(names), arg) != std::end(names);
    };
    if (!is_known())
      throw std::invalid_argument("unknown option " + arg);
    if (i + 1 >= args.size())
      throw std::invalid_argument("option " + arg + " needs a value");
    const std::string &value = args[++i];
    const auto num = [&] { return detail::number(value, arg); };
    if (arg == "-r" || arg == "--radius" || arg == "--telescope-radius") {
      o.telescope_radius_m = num();
    } else if (arg == "-x") {
      o.impact_x_m = num();
    } else if (arg == "-y") {
      o.impact_y_m = num();
    } else if (arg == "-h" || arg == "--h-start") {
      o.h_start_m = num();
    } else if (arg == "--theta" || arg == "--zenith-angle") {
      o.theta_deg = num();
    } else if (arg == "--phi" || arg == "--orientation" || arg == "--azimuth") {
      o.phi_deg = num();
    } else if (arg == "-n" || arg == "--n1" || arg == "--refraction") {
      const double v = num();
      if (v > 1.0 && v < 1.0003)
        o.index_minus_one = v - 1.0;
      else if (v > 0.0 && v < 0.0003)
        o.index_minus_one = v;
      else
        throw std::invalid_argument("invalid refraction index: " + value);
    } else if (arg == "-A" || arg == "--atmosphere") {
      if (value != "const" && value != "-1")
        throw std::invalid_argument("atmosphere tables are not supported: use const or -1");
    } else if (arg == "--index-profile") {
      std::ifstream file(value);
      if (!file)
        throw std::invalid_argument("cannot open index profile " + value);
      o.index_profile = read_index_profile(file, value);
    } else if (arg == "-a" || arg == "--altitude" || arg == "--obslev" ||
               arg == "--observation-level") {
      o.observation_level_m = num();
    } else if (arg == "-b" || arg == "--bs" || arg == "--bunch-size") {
      o.bunch_size = num();
    } else if (arg == "-l" || arg == "--wlmin") {
      o.wavelength_min_nm = num();
    } else if (arg == "-L" || arg == "--wlmax") {
      o.wavelength_max_nm = num();
    } else if (arg == "-q" || arg == "--charge") {
      o.charge = num();
    } else if (arg == "-m" || arg == "--mass") {
      o.mass_gev = num();
    } else if (arg == "-p" || arg == "--momentum") {
      const double p = num();
      if (!(p > 0.0))
        throw std::invalid_argument("momentum must be > 0");
      o.momentum_gev = p;
    } else if (arg == "--beta") {
      beta_given = num();
    } else if (arg == "--primary") {
      o.primary_id = (!value.empty() && std::isdigit(static_cast<unsigned char>(value[0])))
                         ? static_cast<int>(detail::integer(value, arg))
                         : primary_id_from_name(value);
      mass_charge(o.primary_id, o.mass_gev, o.charge);
    } else if (arg == "-N" || arg == "--events") {
      o.events = event_count(detail::integer(value, arg));
    } else if (arg == "--run") {
      o.run_number = detail::integer(value, arg);
    } else if (arg == "--seed") {
      if (!value.empty() && value[0] == '-')
        o.seed = static_cast<std::uint64_t>(detail::integer(value, arg));
      else {
        detail::integer(value, arg);
        o.seed = std::strtoull(value.c_str(), nullptr, 10);
      }
    } else if (arg == "-o" || arg == "--output") {
      parsed.output = value;
    } else if (arg == "--rscat") {
      if (num() != 0.0)
        throw std::invalid_argument("--rscat is not supported");
    } else if (arg == "--nscat") {
      if (detail::integer(value, arg) != 1)
        throw std::invalid_argument("--nscat must be 1 (each event is used once)");
    }
  }
  if (beta_given && o.momentum_gev)
    throw std::invalid_argument("--beta and --momentum are exclusive");
  if (beta_given)
    o.beta = *beta_given;
  if (o.momentum_gev)
    o.beta =
        *o.momentum_gev / std::sqrt(*o.momentum_gev * *o.momentum_gev + o.mass_gev * o.mass_gev);
  return parsed;
}

// Writes the stream of all events; the summary holds rows, weight and the expected weight.
inline StreamSummary write(std::ostream &out, const Options &o) {
  const std::vector<Track> segments = tracks(o);
  RunInfo run;
  run.run_number = o.run_number;
  run.observation_level_m = o.observation_level_m;
  run.wavelength_min_nm = o.wavelength_min_nm;
  run.wavelength_max_nm = o.wavelength_max_nm;
  run.telescopes.push_back({1, {0.0, 0.0, 0.0}, o.telescope_radius_m});
  RunWriter writer(out, "obdeect sources muon-track", run);
  const double energy =
      (o.beta < 1.0 && o.mass_gev > 0.0) ? o.mass_gev / std::sqrt(1.0 - o.beta * o.beta) : 1e-9;
  double expected = 0.0;
  for (const auto &t : segments)
    expected += expected_track_weight(t);
  std::vector<Photon> photons;
  for (std::uint64_t event = 1; event <= o.events; ++event) {
    writer.begin_event(
        {static_cast<std::int64_t>(event), o.primary_id, energy, o.impact_x_m, o.impact_y_m});
    writer.begin_telescope(run.telescopes[0]);
    for (std::size_t k = 0; k < segments.size(); ++k) {
      photons.clear();
      emit_track(segments[k], source_seed(o.seed, k), static_cast<std::uint32_t>(event), 1,
                 photons);
      for (const auto &p : photons)
        writer.push(p);
    }
    writer.end_telescope();
    writer.end_event();
  }
  writer.finish();
  StreamSummary summary = writer.summary();
  summary.expected = expected * static_cast<double>(o.events);
  return summary;
}

} // namespace obdeect::sources::muon_track
