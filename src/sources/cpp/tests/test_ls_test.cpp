#include "check.hpp"

#include "dlamb/programs/ls_test.hpp"
#include "dlamb/programs/octo.hpp"

#include "gstrickd/stream.hpp"

#include <cmath>
#include <cstdint>
#include <sstream>

using namespace dlamb;

namespace {

// ls-test is octo with 4e7 instead of 8e7 photons per position (ls-test.cc:53 and octo.cc:79).
void half_of_octo() {
  const ls_test::Options o;
  octo::Options octo_default;
  octo_default.photons_per_position = 2.0 * o.photons_per_position;
  std::stringstream a(std::ios::in | std::ios::out | std::ios::binary);
  std::stringstream b(std::ios::in | std::ios::out | std::ios::binary);
  const StreamSummary ls = ls_test::write(a, o);
  const StreamSummary oc = octo::write(b, octo_default);
  CHECK_NEAR(ls.expected, 0.5 * oc.expected, 1e-9 * oc.expected);
  CHECK(std::abs(static_cast<double>(ls.photons) - ls.expected) <= 4.0 * std::sqrt(ls.expected));

  gstrickd::Reader reader(a);
  gstrickd::Record record;
  std::uint64_t rows = 0;
  while (reader.next(record))
    if (record.type == gstrickd::RecordType::photons)
      rows += record.chunk.rows();
  CHECK(reader.ended());
  CHECK(rows == ls.photons);
}

} // namespace

int main() {
  return run_checks([] { half_of_octo(); });
}
