#include "check.hpp"

#include "obdeect/stream/stages.hpp"

#include <array>

int main() {
  return run_checks([] {
    // The values are fixed by existing results; a change here changes the random numbers.
    CHECK(obdeect::stream::stage_sources == 0x444C0001U);
    CHECK(obdeect::stream::stage_programs == 0x444C0002U);
    CHECK(obdeect::stream::stage_ray_tracing == 0x444C0003U);
    CHECK(obdeect::stream::stage_atmosphere == 0x444C5546U);
    CHECK(obdeect::stream::stage_ids.size() == 4);
    CHECK(obdeect::stream::detail::all_distinct(obdeect::stream::stage_ids));
    CHECK(!obdeect::stream::detail::all_distinct(std::array<std::uint32_t, 2>{7U, 7U}));
  });
}
