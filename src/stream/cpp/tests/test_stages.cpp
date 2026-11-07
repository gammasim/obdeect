#include "check.hpp"

#include "gstrickd/stages.hpp"

#include <array>

int main() {
  return run_checks([] {
    // The values are fixed by existing results; a change here changes the random numbers.
    CHECK(gstrickd::stage_sources == 0x444C0001U);
    CHECK(gstrickd::stage_programs == 0x444C0002U);
    CHECK(gstrickd::stage_obdeect == 0x444C0003U);
    CHECK(gstrickd::stage_atmosphere == 0x444C5546U);
    CHECK(gstrickd::stage_ids.size() == 4);
    CHECK(gstrickd::detail::all_distinct(gstrickd::stage_ids));
    CHECK(!gstrickd::detail::all_distinct(std::array<std::uint32_t, 2>{7U, 7U}));
  });
}
