#include "check.hpp"

#include "obdeect/sources/version.hpp"

// One version for all components: the sources core reports the package version of the stream
// headers it is built with.
int main() {
  return run_checks([] {
    CHECK(obdeect::sources::version_major == obdeect::stream::version_major);
    CHECK(obdeect::sources::version_minor == obdeect::stream::version_minor);
    CHECK(obdeect::sources::version_patch == obdeect::stream::version_patch);
  });
}
