#include "check.hpp"

#include "obdeect/propagation/streams.hpp"
#include "obdeect/stream/format.hpp"
#include "obdeect/stream/stages.hpp"

int main() {
  return run_checks([] {
    // obdeect propagation is written for photon-stream format major 1 (see streams.hpp).
    CHECK(obdeect::stream::format_major == 1);
    CHECK(obdeect::stream::format_major == obdeect::propagation::expected_gstrickd_format_major);
    // obdeect propagation draws its random numbers from obdeect stream's atmosphere stage id.
    CHECK(obdeect::propagation::stage_atmosphere == obdeect::stream::stage_atmosphere);
  });
}
