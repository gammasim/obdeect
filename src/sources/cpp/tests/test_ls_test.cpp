#include "check.hpp"

#include "dlamb/programs/ls_test.hpp"
#include "dlamb/programs/octo.hpp"
#include "program_stream.hpp"

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

// Emission altitude: every photon of the stream has zem = observation level + z (format.md).
void emission_altitude() {
  const ls_test::Options o;
  std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
  const StreamSummary summary = ls_test::write(stream, o);
  CHECK(summary.photons > 0);
  const StreamContents s = read_stream(stream.str());
  CHECK(!s.rows.empty());
  // ls-test uses the octo observation level (octo.hpp, alt_m).
  const double observation_level_m = octo::Options{}.observation_level_m;
  for (const auto &r : s.rows) {
    CHECK(r.has_zem);
    CHECK_NEAR(r.zem - r.z, observation_level_m, 1e-2);
  }
}

} // namespace

int main() {
  return run_checks([] {
    half_of_octo();
    emission_altitude();
  });
}
