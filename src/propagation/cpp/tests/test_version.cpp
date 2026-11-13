#include "check.hpp"

#include "dlufd/version.hpp"

int main() {
  return run_checks([] {
    CHECK(dlufd::version_major >= 0);
    CHECK(dlufd::version_minor >= 0);
    CHECK(dlufd::version_patch >= 0);
  });
}
