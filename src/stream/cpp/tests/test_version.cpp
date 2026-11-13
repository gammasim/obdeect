#include "check.hpp"

#include "gstrickd/version.hpp"

int main() {
  return run_checks([] {
    CHECK(gstrickd::version_major >= 0);
    CHECK(gstrickd::version_minor >= 0);
    CHECK(gstrickd::version_patch >= 0);
  });
}
