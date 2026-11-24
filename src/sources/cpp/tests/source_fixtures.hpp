#pragma once

// Shared test helpers of the sources tests. Test-only.

#include <cstdio>
#include <exception>

namespace obdeect_test {

// True if action() runs without throwing; a thrown exception is reported and gives false.
template <class Action> bool accepts(Action &&action) {
  try {
    action();
    return true;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "accepts: unexpected exception: %s\n", error.what());
  }
  return false;
}

} // namespace obdeect_test
