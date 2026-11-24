#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string_view>

// Assertions that stay active in Release builds.
inline int check_failures = 0;

#define CHECK(condition)                                                                           \
  do {                                                                                             \
    if (!(condition)) {                                                                            \
      std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #condition);           \
      ++check_failures;                                                                            \
    }                                                                                              \
  } while (false)

#define CHECK_NEAR(value, expected, tolerance)                                                     \
  do {                                                                                             \
    const double check_v = (value), check_e = (expected);                                          \
    if (!(std::abs(check_v - check_e) <= (tolerance))) {                                           \
      std::fprintf(stderr, "%s:%d: CHECK_NEAR failed: %s = %.17g, expected %.17g +- %g\n",         \
                   __FILE__, __LINE__, #value, check_v, check_e, static_cast<double>(tolerance));  \
      ++check_failures;                                                                            \
    }                                                                                              \
  } while (false)

// True if action() throws an Error (std::invalid_argument by default) whose message contains
// `fragment` (any message when empty). Another exception type, or no exception, gives false and
// prints what was thrown, so a test cannot pass on an unrelated error.
template <class Error = std::invalid_argument, class Action>
bool throws_with(Action &&action, std::string_view fragment = {}) {
  try {
    action();
  } catch (const Error &error) {
    if (std::string_view(error.what()).find(fragment) != std::string_view::npos)
      return true;
    std::fprintf(stderr, "throws_with: message \"%s\" lacks \"%.*s\"\n", error.what(),
                 static_cast<int>(fragment.size()), fragment.data());
    return false;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "throws_with: unexpected exception type: %s\n", error.what());
    return false;
  }
  std::fprintf(stderr, "throws_with: no exception\n");
  return false;
}

// Run the tests; an escaping exception counts as a failure.
template <class Tests> int run_checks(Tests tests) noexcept {
  try {
    tests();
  } catch (const std::exception &error) {
    std::fprintf(stderr, "unexpected exception: %s\n", error.what());
    ++check_failures;
  } catch (...) {
    std::fprintf(stderr, "unexpected exception\n");
    ++check_failures;
  }
  if (check_failures)
    std::fprintf(stderr, "%d check(s) failed\n", check_failures);
  return check_failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
