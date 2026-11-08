#include "check.hpp"

#include "dlamb/programs/stream_output.hpp"
#include "dlamb/source.hpp"
#include "dlamb/streams.hpp"

#include "gstrickd/format.hpp"
#include "gstrickd/photons.hpp"
#include "gstrickd/stages.hpp"

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <string>

// dlamb uses gstrickd's stage ids, frame codes, column names and chunk limit directly.
// These checks fail when gstrickd changes a value that dlamb depends on.
int main() {
  return run_checks([] {
    // Format version written and read by dlamb (one named constant, see streams.hpp).
    CHECK(gstrickd::format_major == dlamb::photon_stream_format_major);

    // Stage ids are part of the random-number results: never changed.
    CHECK(dlamb::stage_sources == 0x444C0001U);
    CHECK(dlamb::programs::stage_programs == 0x444C0002U);
    CHECK(gstrickd::stage_sources == dlamb::stage_sources);
    CHECK(gstrickd::stage_programs == dlamb::programs::stage_programs);

    // Frame codes accepted by the Python binding (sky 0, T 2, C 3).
    CHECK(static_cast<std::uint32_t>(gstrickd::Frame::sky) == 0U);
    CHECK(static_cast<std::uint32_t>(gstrickd::Frame::telescope) == 2U);
    CHECK(static_cast<std::uint32_t>(gstrickd::Frame::camera) == 3U);

    // Columns dlamb emits, in gstrickd's names and order.
    const std::string expected[] = {"x", "y", "z", "cx", "cy", "cz", "t", "w", "wavelength"};
    std::size_t k = 0;
    for (std::uint32_t c = 0; c < gstrickd::n_columns; ++c) {
      if ((dlamb::source_columns >> c) & 1U) {
        CHECK(k < std::size(expected));
        if (k < std::size(expected))
          CHECK(std::string(gstrickd::column_name(static_cast<gstrickd::Column>(c))) ==
                expected[k]);
        ++k;
      }
    }
    CHECK(k == std::size(expected));

    CHECK(gstrickd::max_chunk_rows == (1U << 16));
  });
}
