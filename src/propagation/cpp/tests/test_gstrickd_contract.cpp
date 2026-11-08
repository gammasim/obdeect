#include "check.hpp"

#include "dlufd/streams.hpp"
#include "gstrickd/format.hpp"
#include "gstrickd/stages.hpp"

int main() {
  return run_checks([] {
    // dlufd is written for photon-stream format major 1 (see streams.hpp).
    CHECK(gstrickd::format_major == 1);
    CHECK(gstrickd::format_major == dlufd::expected_gstrickd_format_major);
    // dlufd draws its random numbers from gstrickd's atmosphere stage id.
    CHECK(dlufd::stage_atmosphere == gstrickd::stage_atmosphere);
  });
}
