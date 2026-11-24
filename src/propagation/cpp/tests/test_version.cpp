#include "check.hpp"

#include "obdeect/propagation/version.hpp"

// One version for all components: the propagation core reports the package version of the stream
// headers it is built with.
int main() {
  return run_checks([] {
    CHECK(obdeect::propagation::version_major == obdeect::stream::version_major);
    CHECK(obdeect::propagation::version_minor == obdeect::stream::version_minor);
    CHECK(obdeect::propagation::version_patch == obdeect::stream::version_patch);
  });
}
