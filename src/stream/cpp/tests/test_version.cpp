#include "check.hpp"

#include "obdeect/stream/version.hpp"

int main() {
  return run_checks([] {
    CHECK(obdeect::stream::version_major >= 0);
    CHECK(obdeect::stream::version_minor >= 0);
    CHECK(obdeect::stream::version_patch >= 0);
  });
}
