#include "check.hpp"

#include "obdeect/propagation/version.hpp"

int main() {
  return run_checks([] {
    CHECK(obdeect::propagation::version_major >= 0);
    CHECK(obdeect::propagation::version_minor >= 0);
    CHECK(obdeect::propagation::version_patch >= 0);
    CHECK(obdeect::propagation::version_major == obdeect::stream::version_major);
    CHECK(obdeect::propagation::version_minor == obdeect::stream::version_minor);
    CHECK(obdeect::propagation::version_patch == obdeect::stream::version_patch);
  });
}
