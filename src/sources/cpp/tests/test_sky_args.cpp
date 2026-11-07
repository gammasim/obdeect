#include "check.hpp"

#include "../tools/sky_args.hpp"

#include "dlamb/density.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <variant>

using namespace dlamb;
using namespace dlamb::cli;

namespace {

constexpr double deg = pi / 180.0;

// True when fn throws std::invalid_argument (any other exception fails the check).
template <class Fn> bool rejects(Fn fn) {
  try {
    fn();
  } catch (const std::invalid_argument &) {
    return true;
  } catch (...) {
    return false;
  }
  return false;
}

// Event counts: integers in 0..2^32-1; signs, fractions and trailing text are errors.
void counts() {
  CHECK(to_count("0", "-N") == 0U);
  CHECK(to_count("42", "-N") == 42U);
  CHECK(to_count("4294967295", "-N") == std::numeric_limits<std::uint32_t>::max());
  CHECK(rejects([] { (void)to_count("-1", "-N"); }));
  CHECK(rejects([] { (void)to_count("-2147483648", "-N"); }));
  CHECK(rejects([] { (void)to_count("4294967296", "-N"); }));
  CHECK(rejects([] { (void)to_count("1.5", "-N"); }));
  CHECK(rejects([] { (void)to_count("1e3", "-N"); }));
  CHECK(rejects([] { (void)to_count("12x", "-N"); }));
  CHECK(rejects([] { (void)to_count("", "-N"); }));
}

// Fixed wavelength: 200..2000 nm, nothing else.
void fixed_wavelength_range() {
  CHECK(std::get<Fixed>(parse_spectrum("200", "-s")).value == 200.0);
  CHECK(std::get<Fixed>(parse_spectrum("2000", "-s")).value == 2000.0);
  CHECK(rejects([] { (void)parse_spectrum("199.9", "-s"); }));
  CHECK(rejects([] { (void)parse_spectrum("2000.1", "-s"); }));
}

// Cherenkov:from:to is parsed exactly and both ends must lie in 200..2000 nm.
void cherenkov_spectrum() {
  const auto range = std::get<InverseSquare>(parse_spectrum("Cherenkov:300:500", "-s"));
  CHECK(range.low == 300.0 && range.high == 500.0);
  CHECK(std::get<Fixed>(parse_spectrum("cherenkov:300:300", "-s")).value == 300.0);
  CHECK(std::get<InverseSquare>(parse_spectrum("CHERENKOV:200:2000", "-s")).high == 2000.0);
  CHECK(rejects([] { (void)parse_spectrum("Cherenkov:300:500:7", "-s"); }));
  CHECK(rejects([] { (void)parse_spectrum("Cherenkov:300:500xyz", "-s"); }));
  CHECK(rejects([] { (void)parse_spectrum("Cherenkov:300", "-s"); }));
  CHECK(rejects([] { (void)parse_spectrum("Cherenkov:100:500", "-s"); }));
  CHECK(rejects([] { (void)parse_spectrum("Cherenkov:300:2500", "-s"); }));
  CHECK(rejects([] { (void)parse_spectrum("Cherenkov:500:300", "-s"); }));
}

// Pulses: keywords are case-insensitive, parameters must match the arity of the keyword.
void pulse_names_and_arity() {
  const auto simple = std::get<UniformRange>(parse_pulse("SIMPLE:2,5", "-p", false));
  CHECK_NEAR(simple.low, 4.0, 1e-12);
  CHECK_NEAR(simple.high, 6.0, 1e-12);
  const auto simple_lower = std::get<UniformRange>(parse_pulse("simple:2", "-p", false));
  CHECK_NEAR(simple_lower.low, -1.0, 1e-12);
  CHECK_NEAR(simple_lower.high, 1.0, 1e-12);
  CHECK(rejects([] { (void)parse_pulse("Simple:1,2,3", "-p", false); }));
  CHECK(rejects([] { (void)parse_pulse("Simple:", "-p", false); }));

  const auto gauss = std::get<GaussianShape>(parse_pulse("Gauss:2,5", "-p", false));
  CHECK(gauss.sigma == 2.0 && gauss.mean_value == 5.0);
  CHECK(rejects([] { (void)parse_pulse("Gauss:1,2,3", "-p", false); }));
  CHECK(rejects([] { (void)parse_pulse("Exp:1,2,3", "-p", false); }));

  const auto complex = std::get<ComplexPulse>(parse_pulse("Complex:1,2,3", "-p", false));
  CHECK(complex.flat_width == 1.0 && complex.gauss_rms == 2.0 && complex.decay == 3.0);
  CHECK(complex.delay == 0.0);
  const auto complex_delay = std::get<ComplexPulse>(parse_pulse("complex:1,2,3,4", "-p", false));
  CHECK(complex_delay.delay == 4.0);
  CHECK(rejects([] { (void)parse_pulse("Complex:1,2", "-p", false); }));
  CHECK(rejects([] { (void)parse_pulse("Complex:1,2,3,4,5", "-p", false); }));
}

// Angular keywords are case-insensitive; Parallel and Rayleigh were case-sensitive.
void angular_names() {
  CHECK(std::holds_alternative<Isotropic>(parse_angular("Isotropic", "-a")));
  CHECK(std::holds_alternative<Isotropic>(parse_angular("ISOTROPIC", "-a")));
  CHECK(std::holds_alternative<FixedAngle>(parse_angular("PARALLEL", "-a")));
  CHECK(std::holds_alternative<FixedAngle>(parse_angular("parallel", "-a")));
  CHECK(std::holds_alternative<TabulatedAngular>(parse_angular("Rayleigh", "-a")));
  CHECK(std::holds_alternative<TabulatedAngular>(parse_angular("RAYLEIGH", "-a")));
  const auto cone = std::get<FixedAngle>(parse_angular("cone:30", "-a"));
  CHECK_NEAR(cone.u0, 1.0 - std::cos(30.0 * deg), 1e-12);
  CHECK(rejects([] { (void)parse_angular("Cone:30abc", "-a"); }));
}

} // namespace

int main() {
  return run_checks([] {
    counts();
    fixed_wavelength_range();
    cherenkov_spectrum();
    pulse_names_and_arity();
    angular_names();
  });
}
