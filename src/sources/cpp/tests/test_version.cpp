#include "check.hpp"

#include "dlamb/version.hpp"

int main() {
  return run_checks([] {
    CHECK(dlamb::version_major >= 0);
    CHECK(dlamb::version_minor >= 0);
    CHECK(dlamb::version_patch >= 0);
  });
}
