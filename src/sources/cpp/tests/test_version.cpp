#include "check.hpp"

#include "obdeect/sources/version.hpp"

int main() {
  return run_checks([] {
    CHECK(obdeect::sources::version_major >= 0);
    CHECK(obdeect::sources::version_minor >= 0);
    CHECK(obdeect::sources::version_patch >= 0);
  });
}
